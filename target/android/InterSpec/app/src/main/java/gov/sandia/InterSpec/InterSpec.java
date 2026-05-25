/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

package gov.sandia.InterSpec;

import android.app.Activity;
import android.content.res.AssetFileDescriptor;
import android.content.res.AssetManager;
import android.content.res.Configuration;
import android.os.Bundle;

import androidx.appcompat.app.AppCompatActivity;

import android.webkit.*;
import android.net.*;
import android.content.*;
import android.util.*;
import android.view.View;
import android.graphics.Rect;
import android.view.ViewTreeObserver;
import android.os.Message;
import android.os.ParcelFileDescriptor;

import java.io.BufferedInputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileOutputStream;

import android.database.Cursor;
import android.provider.OpenableColumns;
import android.util.Base64;
import android.widget.Toast;

import java.io.InputStream;
import java.security.SecureRandom;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;


/** Native -> Java callback invoked when the C++ side has written a download
 * to a temporary file and we should now show the system Save As dialog. */
interface CallbackFromNativeInterface {
  void callback( String tempFilePath, String displayName );
}

public class InterSpec extends AppCompatActivity
{
  // The Wt server is a process-level singleton -- it binds a TCP port once and stays up for
  // the lifetime of the process, surviving Activity recreations.  Keeping httpPort static is
  // what lets onCreate run on a recreated Activity without trying to start a second server
  // (which would either fail to bind, or succeed and silently leak the old one).
  private static int sHttpPort = 0;
  private static final Object sServerLock = new Object();

  private String interspecID = "";
  private View mDecorView;
  /** File upload callback for platform versions prior to Android 5.0 */
  private ValueCallback<Uri> mUploadMessage;
  /** File upload callback for Android 5.0+ */
  private ValueCallback<Uri[]> mUploadMessageMult;
  private final static int FILECHOOSER_RESULTCODE = 1;
  private final static int OPEN_FILE_CODE = 2;

  private final static int SAVE_AS_CODE = 9821341;

  // When the C++ side finishes writing a download to a temporary file it calls back via JNI;
  // we then launch the system Save As intent and stash the temp file location here until the
  // user picks a destination in onActivityResult.  Initialized to "" so .isEmpty() never NPEs.
  private String mTempFileToSave = "";
  private String mTempFileDisplayName = "";

  private ViewTreeObserver.OnGlobalLayoutListener mLayoutListener;
  private WebView webview;

  // Single-threaded background executor for file I/O that would otherwise ANR on the UI
  // thread: copying picked URIs into the app cache (potentially tens of MB), and writing
  // export bytes back out to the user-chosen save URI.  Serialized so two concurrent
  // exports don't race on the same destination fd.
  private final ExecutorService mIoExecutor = Executors.newSingleThreadExecutor();

  /** Define a class whose member functions we can call from JavaScript. */
  class CallbackFromJavaScriptInterface {
    private InterSpec activity;

    public CallbackFromJavaScriptInterface(InterSpec activity) {
      this.activity = activity;
    }

    /* Function to call when the user wants to open a spectrum file */
    @JavascriptInterface
    public void startBrowseToOpenFile(){
      Log.d("js_inter", "startBrowseToOpenFile()");
      Intent i = new Intent(Intent.ACTION_GET_CONTENT);
      i.addCategory(Intent.CATEGORY_OPENABLE);
      i.setType("*/*");
      activity.startActivityForResult( Intent.createChooser(i,"Select File"),
              OPEN_FILE_CODE );
    }

    // TODO: other places we open up files in InterSpec, define a cooresponding function here to browse for the file and open it; where it makes sense.
  }//class CallbackFromJavaScriptInterface

  public final class WtWebChromeClient extends WebChromeClient 
  {
    @Override
    public boolean onJsAlert(WebView view, String url,  String message,  JsResult result) 
    {
      Log.d("onJsAlert", message);
      return true;
    }

    //openFileChooser for Android 2.2 (API level 8) up to Android 2.3 (API level 10)
	  public void openFileChooser(ValueCallback<Uri> uploadMsg)
	  {
      mUploadMessageMult = null;
	    mUploadMessage = uploadMsg;
	    Intent i = new Intent(Intent.ACTION_GET_CONTENT);
      i.addCategory(Intent.CATEGORY_OPENABLE);
      i.setType("*/*");
      InterSpec.this.startActivityForResult(
	    Intent.createChooser(i,"File Chooser"), FILECHOOSER_RESULTCODE);
    }

	  //openFileChooser for Android 3.0 (API level 11) up to Android 4.0 (API level 15)
    public void openFileChooser(ValueCallback<Uri> uploadMsg, String acceptType)
    {
      Log.d("onShowFileChooser", "In chooser for Android 3.0 through 4.0");
      openFileChooser(uploadMsg);
    }

    //openFileChooser for Android 4.1 (API level 16) up to Android 4.3 (API level 18)
    public void openFileChooser(ValueCallback<Uri> uploadMsg, String acceptType, String capture) 
	  {
      Log.d("onShowFileChooser", "In chooser for Android 4.1 through 4.3");
      openFileChooser(uploadMsg);
    }

    //Android 4.4 - out of luck :(

    //Android 5.0 (API level 21) and above
    public boolean onShowFileChooser( WebView webView, ValueCallback<Uri[]> filePathCallback, WebChromeClient.FileChooserParams fileChooserParams )
    {
      Log.d("onShowFileChooser", "In chooser for Android >5.0");
      mUploadMessageMult = filePathCallback;
      mUploadMessage = null;
      Intent i = new Intent(Intent.ACTION_GET_CONTENT);
      i.addCategory(Intent.CATEGORY_OPENABLE);
      i.setType("*/*");
      InterSpec.this.startActivityForResult(
              Intent.createChooser(i,"File Chooser"), FILECHOOSER_RESULTCODE);

      return true;
    }

    /** Intercept target='_blank' / window.open() navigations so we can route Wt download links
     *  through the native save path.  Wt sets target='_blank' on its WLink-as-anchor downloads,
     *  which would otherwise produce a do-nothing onCreateWindow with no URL exposed.  Trick:
     *  create a throwaway WebView, attach a WebViewClient whose shouldOverrideUrlLoading captures
     *  the URL, and hand the new view back to the framework so the load actually starts.  We
     *  then route that captured URL through nativeDownloadUrl() and abort the load. */
    @Override
    public boolean onCreateWindow( WebView view, boolean isDialog, boolean isUserGesture, Message resultMsg )
    {
      Log.d("onCreateWindow", "Intercepting target=_blank navigation");

      final WebView captureView = new WebView( view.getContext() );
      captureView.setWebViewClient( new WebViewClient() {
        @Override
        public boolean shouldOverrideUrlLoading( WebView v, String url )
        {
          Log.d("onCreateWindow", "captured url=" + url);
          handlePossibleDownloadUrl( url );
          v.destroy();
          return true;
        }
      } );

      final WebView.WebViewTransport transport = (WebView.WebViewTransport) resultMsg.obj;
      transport.setWebView( captureView );
      resultMsg.sendToTarget();
      return true;
    }

    @Override
    public boolean onConsoleMessage (ConsoleMessage consoleMessage)
    {
      Log.d("onConsoleMessage", consoleMessage.message());
      return true;
    }
  }//public final class WtWebChromeClient extends WebChromeClient 


  /** Route a URL captured from the WebView (target=_blank or DownloadListener) through the
   *  native download-to-save path.  Only localhost URLs are treated as downloads; anything else
   *  is logged and ignored.  External links could be opened in the system browser in a future
   *  enhancement (compare the wxWidgets target). */
  private void handlePossibleDownloadUrl( String url )
  {
    if( url == null || url.isEmpty() )
      return;

    final Uri u;
    try {
      u = Uri.parse( url );
    } catch( Exception e ) {
      Log.w("download", "Could not parse url='" + url + "': " + e);
      return;
    }

    final String host = (u != null) ? u.getHost() : null;
    if( "127.0.0.1".equals(host) || "localhost".equals(host) )
    {
      Log.d("download", "Routing localhost url through native save: " + url);
      nativeDownloadUrl( url );
    }
    else
    {
      Log.w("download", "Ignoring non-localhost download url='" + url + "' (host='" + host + "')");
    }
  }

  public String getDisplayFileName( Uri result )
  {
    String name = null;
    Cursor cursor = InterSpec.this.getContentResolver().query(result, null, null, null, null, null);
    try 
    {
	  if( cursor != null && cursor.moveToFirst() )
        name = cursor.getString(cursor.getColumnIndexOrThrow(OpenableColumns.DISPLAY_NAME));
    }finally
    {
	  if( cursor != null )
        cursor.close();
    }
	
    return name;
  }//public String getDisplayFileName( Uri result ) 

  /** Copies a content:// URI into the app's cache directory and returns the on-disk path.
   *  IO-heavy and must be called from a worker thread (use mIoExecutor); throws on UI thread
   *  for large files. */
  private String copyUriToTmpDir( Uri result, final String displayName )
  {
    try
    {
      Log.d("copyUriToTmpDir", "About to open ParcelFileDescriptor for " + result );
      final ParcelFileDescriptor parcelFD = getContentResolver().openFileDescriptor(result, "r");
      if( parcelFD == null )
      {
        Log.w("copyUriToTmpDir", "openFileDescriptor returned null for " + result);
        return null;
      }

      final File tmpDir = getCacheDir();
      final File outputFile = File.createTempFile("spectrum", ".tmp", tmpDir);
      outputFile.deleteOnExit();

      try( InputStream input = new FileInputStream( parcelFD.getFileDescriptor() );
           FileOutputStream output = new FileOutputStream( outputFile ) )
      {
        final byte[] buffer = new byte[64 * 1024];
        int length;
        while( (length = input.read(buffer)) > 0 )
          output.write(buffer, 0, length);
        output.flush();
      }
      finally
      {
        parcelFD.close();
      }

      Log.d("copyUriToTmpDir", "Done copying " + result + " -> " + outputFile.getAbsolutePath() );
      return outputFile.getAbsolutePath();
    }
    catch( IOException e )
    {
      Log.w("copyUriToTmpDir", "IOException copying " + result, e);
      runOnUiThread( new Runnable() { public void run() {
        Toast.makeText(getApplicationContext(), "Error copying file into InterSpec.", Toast.LENGTH_SHORT).show();
      } } );
      return null;
    }
    catch( NullPointerException e )
    {
      Log.w("copyUriToTmpDir", "NPE copying " + result, e);
      runOnUiThread( new Runnable() { public void run() {
        Toast.makeText(getApplicationContext(), "Error (2) copying file into InterSpec.", Toast.LENGTH_SHORT).show();
      } } );
      return null;
    }
  }

  /** Dispatches the URI copy + InterSpec open onto the IO executor.  Safe to call from the
   *  UI thread; large spectrum files (tens of MB) won't ANR. */
  public void openFileInInterSpec( final Uri result )
  {
    Log.d("openFileInInterSpec", "Got URI and interspecID='" + interspecID + "', httpPort=" + sHttpPort );

    if( result == null )
      return;

    final String scheme = (result.getScheme() != null) ? result.getScheme().toLowerCase() : "";
    if( "interspec".equals(scheme) || "raddata".equals(scheme) )
    {
      Log.d("openFileInInterSpec", "Scheme=" + scheme + ", url=" + result );
      openAppUrl( interspecID, result.toString() );
      return;
    }

    mIoExecutor.execute( new Runnable() { public void run() {
      String displayName = getDisplayFileName( result );
      if( displayName == null )
        displayName = "unnamedfile";

      final String pathname = copyUriToTmpDir( result, displayName );
      if( pathname == null )
        return;

      Log.d("openFileInInterSpec", "Will send to InterSpec: " + pathname);
      final int status = openFile( interspecID, pathname );
      if( status <= 0 )
        Log.w("openFileInInterSpec", "openFile failed with status=" + status + " for " + pathname);
      // The temp file is deleteOnExit-marked, so no explicit cleanup here.
    } } );
  }

  @Override
  protected void onActivityResult( int requestCode, int resultCode, Intent intent )
  {
    super.onActivityResult(requestCode, resultCode, intent);

    if( requestCode == SAVE_AS_CODE )
    {
      final String tempPath = mTempFileToSave;
      mTempFileToSave = "";
      mTempFileDisplayName = "";

      if( resultCode != Activity.RESULT_OK )
      {
        // User cancelled.  Clean up the temp file C++ wrote for them.
        if( !tempPath.isEmpty() )
        {
          final File file = new File( tempPath );
          Log.d("SAVE_AS", "User cancelled; deleted temp=" + file.delete() + " path=" + tempPath );
        }
        return;
      }

      if( tempPath.isEmpty() )
      {
        Toast.makeText( getApplicationContext(),
                        "No output file to copy to destination - sorry!",
                        Toast.LENGTH_SHORT ).show();
        return;
      }

      final Uri destUri = intent.getData();
      if( destUri == null )
      {
        Log.w("SAVE_AS", "RESULT_OK but intent.getData() was null");
        return;
      }

      // Copy temp -> user-chosen URI on the IO executor; toast back from the UI thread.
      // Large spectrum/N42 exports can be tens of MB so doing this synchronously on the UI
      // thread would ANR.
      mIoExecutor.execute( new Runnable() { public void run() {
        final Context context = getApplicationContext();
        long totalBytes = 0;
        boolean ok = false;

        try( OutputStream output = context.getContentResolver().openOutputStream(destUri);
             FileInputStream is = new FileInputStream(tempPath) )
        {
          if( output == null )
            throw new IOException("openOutputStream returned null for " + destUri);

          final byte[] buffer = new byte[64 * 1024];
          int length;
          while( (length = is.read(buffer)) > 0 )
          {
            output.write(buffer, 0, length);
            totalBytes += length;
          }
          output.flush();
          ok = true;
          Log.d("SAVE_AS", "Wrote " + totalBytes + " bytes from " + tempPath + " -> " + destUri);
        }
        catch( IOException e )
        {
          Log.w("SAVE_AS", "Failed writing " + destUri, e);
        }
        finally
        {
          final File tmp = new File(tempPath);
          if( !tmp.delete() )
            Log.w("SAVE_AS", "Failed to delete temp file " + tempPath );
        }

        final boolean success = ok;
        runOnUiThread( new Runnable() { public void run() {
          Toast.makeText( context,
                          success ? "Done writing file." : "Error writing output file.",
                          Toast.LENGTH_SHORT ).show();
        } } );
      } } );

      return;
    }//if( this is result of save as dialog )

    if( requestCode == FILECHOOSER_RESULTCODE )
    {
      // The single-URI callback (pre-API 21) wants exactly one Uri.
      if( mUploadMessage != null )
      {
        final Uri result = (intent == null || resultCode != RESULT_OK) ? null : intent.getData();
        Log.d("onActivityRe", "single-upload result=" + result);
        openFileInInterSpec( result );
        mUploadMessage.onReceiveValue( result );
        mUploadMessage = null;
      }
      // The multi-URI callback (API 21+) accepts any number of Uris.  Multi-select picks come
      // back in intent.getClipData(); a single-file pick comes back in intent.getData().  The
      // original code only ever looked at getData(), so even when the user selected three
      // files in the system picker, only the first one made it through.
      else if( mUploadMessageMult != null )
      {
        Uri[] dataUris = null;
        if( intent != null && resultCode == RESULT_OK )
        {
          final List<Uri> uris = new ArrayList<>();
          final ClipData clip = intent.getClipData();
          if( clip != null )
          {
            for( int i = 0; i < clip.getItemCount(); i++ )
            {
              final Uri u = clip.getItemAt(i).getUri();
              if( u != null )
                uris.add(u);
            }
          }
          else if( intent.getData() != null )
          {
            uris.add( intent.getData() );
          }
          if( !uris.isEmpty() )
            dataUris = uris.toArray(new Uri[0]);
        }

        Log.d("onActivityRe", "multi-upload count=" + (dataUris == null ? 0 : dataUris.length));
        if( dataUris != null )
        {
          for( Uri u : dataUris )
            openFileInInterSpec(u);
        }
        mUploadMessageMult.onReceiveValue( dataUris );
        mUploadMessageMult = null;
      }
      else
      {
        Log.w("onActivityRe", "FILECHOOSER_RESULTCODE with neither upload callback set");
      }
    }//if( requestCode == FILECHOOSER_RESULTCODE )

    if (requestCode == OPEN_FILE_CODE) {
      Log.d("onActivityRe", "Picked file based on starting from JS");

      Uri result = intent == null || resultCode != RESULT_OK ? null : intent.getData();

      openFileInInterSpec(result);
    }
  }//protected void onActivityResult(...)

  @Override
  public void onConfigurationChanged(Configuration newConfig) {
    super.onConfigurationChanged(newConfig);
    // The configChanges= attribute in AndroidManifest.xml lists everything we want to handle
    // without recreating the Activity; the WebView itself responds to layout changes via
    // the OnGlobalLayoutListener registered in onStart, so there's nothing to do here.
  }

  @Override
  public void onCreate( Bundle savedInstanceState )
  {
    // super.onCreate must run before anything else so AppCompatActivity wires up properly;
    // if any conditional below threw before we called super, the framework would raise
    // SuperNotCalledException and the app would crash on every config change.
    super.onCreate(savedInstanceState);

    Log.d("onCreate", "Starting");

    Intent intent = getIntent();
    String action = intent.getAction();
    String type = intent.getType();
    Uri data = intent.getData();

    if( data != null )
      Log.d("onCreate", "data: " + data.toString() );
    else
      Log.d("onCreate", "data: null" );

    // SecureRandom + 128 bits of entropy: the Wt server binds to 127.0.0.1, which on Android
    // is shared between all installed apps, so any app on the device can poke at the port.
    // The original java.util.Random + nextInt(9999999) gave ~10^7 guesses and a predictable
    // time-based seed -- trivially brute-forceable.  Base64 URL-safe encoding stays inside
    // the character set Wt's session-token validator accepts.
    final SecureRandom secureRandom = new SecureRandom();
    final byte[] tokenBytes = new byte[16];
    secureRandom.nextBytes(tokenBytes);
    interspecID = "session" + Base64.encodeToString(tokenBytes,
                                                    Base64.URL_SAFE | Base64.NO_WRAP | Base64.NO_PADDING);

    addPrimarySessionToken( interspecID );

    // Re-register the JNI -> Java save-file callback on every onCreate so it captures *this*
    // Activity instance.  On Activity recreation the old anonymous class instance closes over
    // the previous (now destroyed) Activity and runOnUiThread / startActivityForResult on it
    // are no-ops or crashes.  setFileSaveCallback replaces the previous global ref C++-side.
    setFileSaveCallback( new CallbackFromNativeInterface() {
      public void callback( final String tempFilePath, final String displayName ) {
        Log.d("fileSaveCallback", "tempFilePath='" + tempFilePath + "', displayName='" + displayName + "'");

        if( tempFilePath == null || tempFilePath.isEmpty() ) {
          runOnUiThread(new Runnable() {
            public void run() {
              Toast.makeText( getApplicationContext(),
                              "Internal error: native download produced no file.",
                              Toast.LENGTH_SHORT ).show();
            }
          });
          return;
        }

        runOnUiThread(new Runnable() {
          public void run() {
            mTempFileToSave = tempFilePath;
            mTempFileDisplayName = (displayName != null) ? displayName : "download";

            Intent intent = new Intent(Intent.ACTION_CREATE_DOCUMENT);
            intent.addCategory(Intent.CATEGORY_OPENABLE);
            intent.setType("*/*");
            intent.putExtra(Intent.EXTRA_TITLE, mTempFileDisplayName);
            startActivityForResult(intent, SAVE_AS_CODE);
          }
        });
      }
    });

    // Start the Wt server exactly once per process.  sHttpPort survives Activity recreation
    // (config change, theme switch, etc.) so we don't try to bind a new port -- which would
    // fail if the previous instance's port is still held, or worse, succeed and leak it.
    synchronized( sServerLock )
    {
      if( sHttpPort == 0 )
      {
        Log.i("onCreate", "First-time server start");
        initNative();
        setTempDir( getCacheDir().getPath() );
        sHttpPort = startWt(this);
      }
      else
      {
        Log.i("onCreate", "Server already running on port " + sHttpPort);
      }
    }

    // No-title / no-action-bar is set declaratively via Theme.AppCompat.DayNight.NoActionBar
    // in res/values/styles.xml (AppTheme), so no programmatic requestWindowFeature is needed.

    {
      setContentView(R.layout.main);

      webview = (WebView)findViewById(R.id.webview);
      WebSettings settings = webview.getSettings();
      // LOAD_DEFAULT lets the WebView cache Wt's static assets (CSS, JS, fonts), which on a
      // cold start is the difference between an instant render and a multi-second hourglass.
      // Wt versions assets with their fingerprints so cache invalidation happens naturally.
      settings.setCacheMode(WebSettings.LOAD_DEFAULT);
      // Must be true so that target='_blank' navigations (which Wt sets on its download
      // anchors) reach WebChromeClient.onCreateWindow, where we capture the URL and route it
      // through nativeDownloadUrl().  If false, the WebView silently drops them.
      settings.setSupportMultipleWindows(true);
      settings.setJavaScriptCanOpenWindowsAutomatically(true);
      settings.setJavaScriptEnabled(true);
      //settings.setLoadWithOverviewMode(true);
      //settings.setUseWideViewPort(true);
      //settings.setSupportZoom(true);
      //settings.setGeolocationEnabled(false);
      //settings.setBuiltInZoomControls(false);
      //settings.setLayoutAlgorithm(WebSettings.LayoutAlgorithm.SINGLE_COLUMN);
      //settings.setAllowFileAccess(true);  // TODO: explore being able to save files directly from the WebView - maybe this is only accesing existing files?
      settings.setDomStorageEnabled(true);
      webview.setScrollBarStyle(WebView.SCROLLBARS_OUTSIDE_OVERLAY);
      webview.setScrollbarFadingEnabled(true);
      webview.setLayerType(View.LAYER_TYPE_HARDWARE, null);
      webview.setVerticalScrollBarEnabled(false);
      webview.setHorizontalScrollBarEnabled(false);
      webview.setScrollContainer(false);

      CallbackFromJavaScriptInterface jsInterface = new CallbackFromJavaScriptInterface(this);
      webview.addJavascriptInterface(jsInterface, "interspecJava");

      webview.setWebChromeClient( new WtWebChromeClient() );
      // Belt-and-suspenders download intercept: if Wt ever serves a resource that hits the WebView's
      // download manager directly (Content-Disposition: attachment on a same-window navigation),
      // capture it here and route through the native save path.  The primary path for Wt downloads
      // is target='_blank' anchors, which arrive via WebChromeClient.onCreateWindow above.
      webview.setDownloadListener( new DownloadListener() {
        @Override
        public void onDownloadStart( String url, String userAgent, String contentDisposition,
                                     String mimetype, long contentLength )
        {
          Log.d("onDownloadStart", "url=" + url + ", length=" + contentLength
                                   + ", contentDisposition='" + contentDisposition + "'");
          handlePossibleDownloadUrl( url );
        }
      } );

      webview.setFocusable(true);
      webview.setFocusableInTouchMode(true);
      webview.requestFocus(View.FOCUS_DOWN);

      WebViewClient ourWebClient = new WebViewClient(){
        @Override
        public boolean shouldOverrideUrlLoading( WebView view, String url )
        {
          // Wt is an all-Ajax framework: legitimate page state never arrives via a top-level
          // navigation in our WebView -- only the initial webview.loadUrl() does that, and
          // shouldOverrideUrlLoading is not consulted for it.  Any other URL load here is a
          // Wt internal anchor or form submit that would replace the current page and
          // re-trigger saved-state restoration (which can fire long-running fits, etc.).
          // Returning true blocks that.  Downloads still work via onCreateWindow (target=_blank)
          // and the DownloadListener.
          Log.d("shouldOverrideUrlLoad", "blocking top-level nav to url=" + url);
          return true;
        }
      };
      webview.setWebViewClient( ourWebClient );



      if( Intent.ACTION_VIEW.equals(action)
              || Intent.ACTION_EDIT.equals(action)
              || Intent.ACTION_PICK.equals(action)
              || Intent.ACTION_DEFAULT.equals(action)
              //|| Intent.ACTION_BROWSABLE.equals(action)
       )
      {
        Log.d("onCreate", "ACTION_VIEW: Will set initial file to load at start");
        Uri fileUri = (Uri) intent.getData();

        if( fileUri != null )
        {
          String scheme = fileUri.getScheme();
          Log.d("onCreate", "scheme: '" + scheme + "'" );

          if( scheme.equalsIgnoreCase("raddata") || scheme.equalsIgnoreCase("interspec") )
          {
            int status = setInitialFileToLoad( interspecID, fileUri.toString() );
            Log.d("onCreate", "ACTION_VIEW, Will open URI on load; status=" + status );
          }else if( ContentResolver.SCHEME_CONTENT.equals(scheme) )
          {
            // handle as content uri
            Log.d("onCreate", "ACTION_VIEW, fileUri is a SCHEME_CONTENT" );

            String displayName = getDisplayFileName( fileUri );
            if( displayName == null )
              displayName = "unamedfile";
            String pathname = copyUriToTmpDir( fileUri, displayName );
            if( pathname != null )
            {
              int status = setInitialFileToLoad( interspecID, pathname );
              Log.d("onCreate", "ACTION_VIEW, Will open file '" + pathname + "' on load; status=" + status );
              //We should delete the file at some point...
            }
          }else
          {
            // handle as file uri
            int status = setInitialFileToLoad( interspecID, fileUri.getPath() );
            Log.d("onCreate", "ACTION_VIEW, fileUri is NOT a SCHEME_CONTENT, but a path=" + fileUri.getPath()  + "; will open with status=" + status );
          }
        }else
        {
          Log.d("onCreate", "ACTION_VIEW, fileUri is null");
        }
      }


      else if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null )
      {
        Log.d("onCreate", "ACTION_SEND_MULTIPLE");
        final ArrayList<Uri> fileUris = intent.getParcelableArrayListExtra(Intent.EXTRA_STREAM);
        if( fileUris != null )
        {
          for( Uri u : fileUris )
            openFileInInterSpec(u);
        }
      }else
      {
        Log.d("onCreate", "Another Action: " + action );
      }

      Log.d("onCreate", "Will load 'http://127.0.0.1:" + sHttpPort + "/?apptoken=" + interspecID + "'");
      webview.loadUrl("http://127.0.0.1:" + sHttpPort + "/?apptoken=" + interspecID );

      // WebView.setWebContentsDebuggingEnabled(true);
    }

    Log.d("onCreate", "done starting server ish");

    mDecorView = getWindow().getDecorView();

    final View mActivityRootView = findViewById(android.R.id.content);

    mLayoutListener = new ViewTreeObserver.OnGlobalLayoutListener() {

      /// When the soft-keyboard shows, we'll query the position of the input element
      /// asynchronously using JavaScript, so when we get a result, we will store it
      /// in this next variable, then request a layout, and then move things when
      /// that translate the WebView when that happens (it looks like we can only translate
      /// the WebView after a layout cycle completes).
      private double mInputYFromJS = 0.0;

      // We'll track the keyboard showing, but we dont actually use this right now
      private boolean mShowingKeyboard = false;
      /// I cant get Android to move the WebView so that the input is visible, and we dont want to
      /// resize the WebView when the soft keyboard shows, so we will watch for layout changes, and
      /// detect the keyboard there.
      /// For Android 11 (API level 30), and newer, there is a WindowInsets API that would be better
      /// to use, but this would leave out too many users right now.
      @Override
      public void onGlobalLayout() {
        Rect r = new Rect();
        mActivityRootView.getWindowVisibleDisplayFrame(r);

        int activityH = mActivityRootView.getHeight(); // This is app area, minus top and bottom OS bars.  Not affected if keyboard is shown.
        int visibleH = (r.bottom - r.top); // This is app area, accounting for keyboard being shown (so height doesnt include keyboard height, just the visible part of the app)
        int heightDiff = activityH - visibleH;

        int screenHeight = mActivityRootView.getRootView().getHeight(); // Includes the OS top and bottom bars
        Context context = getApplicationContext();
        //Toast.makeText(context, String.format("screenHeight: %d, activityHeight: %d, visibleHeight: %d, heightDiff: %d", screenHeight, activityH, visibleH, heightDiff), Toast.LENGTH_SHORT ).show();

        // heightDiff is zero, if the keyboard is not showing, but at least a few hundred px if it is showing
        if (heightDiff > 100) {
          mShowingKeyboard = true;

          // We may be here either when the keyboard is first showing, OR as a result of us
          // requesting a layout below when we get the input element position from the WebView via JS
          if( mInputYFromJS > 0.0 ) {
            // If mInputYFromJS is greater than zero, then we are here as a result of getting the
            // element position from the JS, and we can now slide the WebView up the right amount

            // We'll assume element is 20px high, and add another 20px padding
            double posToMakeVisibleDbl = r.top + mInputYFromJS + 20.0 + 20.0;
            int posToMakeVisible = (int) Math.round(posToMakeVisibleDbl);

            //On phone when keyboard shows: screenHeight: 1280, activityH: 1136, visibleH: 634, posToMakeVisible: 743
            //Toast.makeText(context, String.format("screenHeight: %d, activityH: %d, visibleH: %d, posToMakeVisible: %d", screenHeight, activityH, visibleH, posToMakeVisible), Toast.LENGTH_SHORT ).show();

            if( posToMakeVisible > visibleH )
            {
              int delta = visibleH - posToMakeVisible;
              webview.setTranslationY( delta );
            }else
            {
              //Toast.makeText(context, String.format("Setting translation to zero, posToMakeVisible=%d, visibleH=%d", posToMakeVisible, visibleH), Toast.LENGTH_SHORT).show();
              webview.setTranslationY( 0 );
            }
          } else {
            // If mInputYFromJS is zero, then the keyboard just appeared, and we dont have the input
            // element position, and we need to get it.

            // Ask JS for the active input's screen-pixel Y by multiplying its CSS-pixel
            // bounding-rect.top by window.devicePixelRatio in the WebView itself, so we don't
            // have to guess whether to apply the system DisplayMetrics density on the Java
            // side.  The old code used TypedValue.applyDimension(COMPLEX_UNIT_DIP, ...), which
            // only matches when the WebView's devicePixelRatio equals the system density --
            // not true on tablets that set their own viewport meta or initial scale.
            webview.evaluateJavascript(
                    "(function() { " +
                            "  const el = document.activeElement; " +
                            "  if( !el || (el.tagName != 'INPUT') ) return '-1';" +
                            "  const rect = el.getBoundingClientRect(); " +
                            "  return '' + (rect.top * (window.devicePixelRatio || 1));" +
                            "})()",
                    value -> {
                      try {
                        double topPx = Double.parseDouble(value.replaceAll("[\"']", "")); //value is something like "\"519.667\"", so we need to get rid of leading/trailing quote

                        if( topPx < 0.0 ) {
                          mInputYFromJS = topPx; // There was an error.
                        }else{
                          mInputYFromJS = topPx; // Already in device pixels (CSS px * dpr).
                        }

                        // We cant access the WebView from here, and even if we evaluate on the GUI
                        //  thread, thats not good enough to alter the layout, instead we have to
                        //  request the layout to be computed, which will call the onGlobalLayout()
                        //  where we can then set the Y translation
                        webview.requestLayout();
                      } catch (NumberFormatException e) {
                        Log.e("WebView", "Error parsing position", e);
                        Toast.makeText(context, "Error parsing result from webview, returned: '" + value + "', errmsg=" + e.toString(), Toast.LENGTH_SHORT).show();
                      }
                    }); //webview.evaluateJavascript(...)
            webview.setTranslationY( 0 ); // JIC
            mInputYFromJS = 0.0; //JIC
          }//if( mInputYFromJS > 0.0 ) / else
        } else {
          // Reset translation back to zero when keyboard is hidden
          webview.setTranslationY( 0 );
          mInputYFromJS = 0.0;
          mShowingKeyboard = false;
        }
      }
    };

  }//public void onCreate(Bundle savedInstanceState)

  @Override
  protected void onStart() {
    super.onStart();
    mDecorView.getViewTreeObserver().addOnGlobalLayoutListener(mLayoutListener);
  }

  @Override
  protected void onStop() {
    super.onStop();
    mDecorView.getViewTreeObserver().removeOnGlobalLayoutListener(mLayoutListener);
  }

  @Override
  protected void onPause()
  {
    super.onPause();
    Log.d("onPause", "onPause was called");
  }//void onPause()


  @Override
  protected void onResume()
  {
    super.onResume();

    Log.d("onResume", "onResume was called");
  }//void onResume()


  @Override
  protected void onNewIntent( Intent intent )
  {
    Log.d("onNewIntent", "onNewIntent" );

    //This function is called when InterSpec is already running (but in the background)
    //  and another app like Google Drive requests it to open up a file.	  
	super.onNewIntent( intent );
    setIntent( intent );
    String action = intent.getAction();
    String type = intent.getType();

    if( Intent.ACTION_VIEW.equals(action)
            || Intent.ACTION_EDIT.equals(action)
            || Intent.ACTION_PICK.equals(action)
            || Intent.ACTION_DEFAULT.equals(action)
    //|| Intent.ACTION_BROWSABLE.equals(action)
      )
    {
      Uri fileUri = (Uri) intent.getData();

      openFileInInterSpec( fileUri );
    }else if( Intent.ACTION_SEND.equals(action) )
    {
      Uri fileUri = (Uri) intent.getParcelableExtra(Intent.EXTRA_STREAM);
      if( fileUri != null )
        openFileInInterSpec( fileUri );
    }else if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null )
    {
      final ArrayList<Uri> fileUris = intent.getParcelableArrayListExtra(Intent.EXTRA_STREAM);
      if( fileUris != null )
      {
        for( Uri u : fileUris )
          openFileInInterSpec(u);
      }
    } else
    //Intent.ACTION_QUICK_VIEW.equals(action)
    //Intent.ACTION_QUICK_VIEW.equals(action)
    {
        // Handle other intents, such as being started from the home screen
    }//if( decide what actiuon is wanted )
  }//void onNewIntent (Intent intent)
  

  @Override
  protected void onDestroy()
  {
    // Tell the Wt server to drop our session token; otherwise an Activity recreation would
    // register a new token without ever releasing the old one (a server-side leak that adds
    // up over many config changes in a long-running session).
    if( interspecID != null && !interspecID.isEmpty() )
    {
      try {
        removeSessionToken( interspecID );
      } catch( Throwable t ) {
        Log.w("onDestroy", "removeSessionToken failed for " + interspecID, t);
      }
    }
    super.onDestroy();
  }
  
  /*
  @Override
  public void onWindowFocusChanged( boolean hasFocus )
  {
    super.onWindowFocusChanged( hasFocus );

    mHideHandler.removeMessages(0);
    if( hasFocus )
    {
      mHideHandler.sendEmptyMessageDelayed(0, 300);
    }
  }//public void onWindowFocusChanged(boolean hasFocus) 
  */

  /*
  private void hideSystemUI() 
  {
    mDecorView.setSystemUiVisibility(
                  View.SYSTEM_UI_FLAG_LAYOUT_HIDE_NAVIGATION
                  | View.SYSTEM_UI_FLAG_LAYOUT_FULLSCREEN
                  | View.SYSTEM_UI_FLAG_HIDE_NAVIGATION
                  | View.SYSTEM_UI_FLAG_FULLSCREEN
                  | View.SYSTEM_UI_FLAG_LOW_PROFILE
                  | View.SYSTEM_UI_FLAG_IMMERSIVE_STICKY);

	//View.SYSTEM_UI_FLAG_IMMERSIVE View.SYSTEM_UI_FLAG_IMMERSIVE_STICKY
  }//private void hideSystemUI() 
*/

  /*
  private final Handler mHideHandler = new Handler()
  {
    @Override
    public void handleMessage(Message msg) 
	{
      hideSystemUI();
    }
  };
  */

	  
  
/*  
  public void onSaveInstanceState(Bundle savedInstanceState) 
  {
    super.onSaveInstanceState(savedInstanceState);
    
    savedInstanceState.putInt("httpPort", httpPort);
    
//    savedInstanceState.putBoolean("MyBoolean", true);
//    savedInstanceState.putDouble("myDouble", 1.9);
//    savedInstanceState.putInt("MyInt", 1);
//    savedInstanceState.putString("MyString", "Welcome back to Android");
  }
  
  @Override
  public void onRestoreInstanceState(Bundle savedInstanceState) 
  {
    super.onRestoreInstanceState(savedInstanceState);
    
    if( !savedInstanceState.containsKey("httpPort") )
    {
      return;
    }
    
//    boolean myBoolean = savedInstanceState.getBoolean("MyBoolean");
//    httpPort = savedInstanceState.getInt("httpPort");
//    String myString = savedInstanceState.getString("MyString");
  }
*/

  static
  {
    // Load the 'native-lib' library on application startup.
    System.loadLibrary("InterSpecAppLib" );
  }

  public static int startWt( Activity activity )
  {
    Log.d("WtAndroid::startWt", "in startWt ...");

    String wtAssetsDir = activity.getFilesDir().getAbsolutePath() + "/wt-assets/";
    try
    {
      copyWtAssets(wtAssetsDir, activity.getAssets());
    } catch (IOException e)
    {
      Log.d("startWt", "Caught exception extracting interspec-assets.zip, stacktrace:" );
      e.printStackTrace();
    }

    //Note: as of 20190311, the getExternalFilesDir() or getFilesDir() calls have not been tested.
    String userDataDir = activity.getExternalFilesDir(null).getAbsolutePath();
    if( userDataDir.isEmpty() )
      userDataDir = activity.getFilesDir().getAbsolutePath();

    String xml_config_path = wtAssetsDir + "/data/config/wt_config_android.xml";

    // TODO: need to change wt_error_log.log to go into cache directory

    Log.d("WtAndroid::startWt", "About to call into native");
    int httpPort = startServingInterSpec( "InterSpec", userDataDir, wtAssetsDir, xml_config_path );

    if( httpPort <= 0 )
    {
      Log.d("WtAndroid::startWt", "Failed to  " + httpPort);
      // TODO: display flat HTML page with error
    }

    Log.d("WtAndroid::startWt", "Started wt application on http-port " + httpPort);

    return httpPort;
  }//public static int startWt(...)


  private static void copyWtAssets(String wtAssetsDir, AssetManager am) throws IOException, SecurityException
  {
    //We need to know if the assets have been updated since we last ran things to see if we need to
    // re-extract things (takes ~6 seconds to extract).  So as a hack for doing this, we will see
    // if interspect-assets.zip has changed size since last time, and only extract it if it has.
    // I know this isnt quite right, but close enough for now.

    final File parent_dir = new File(wtAssetsDir);
    final File sizeMarker = new File(wtAssetsDir + "/asset_size.txt");

    // Bundled-asset size is obtained from the AssetFileDescriptor in O(1) when the zip is
    // stored uncompressed in the APK (see noCompress 'zip' in app/build.gradle).  Older builds
    // had it compressed -- openFd then throws -- so fall back to a stream-and-skip pass that
    // still beats the original code's "open the zip twice and count" by 2x.
    long bundledZipLen = -1;
    try( AssetFileDescriptor afd = am.openFd("interspec-assets.zip") )
    {
      bundledZipLen = afd.getLength();
    }
    catch( IOException e )
    {
      Log.d("copyWtAssets", "openFd failed (zip is APK-compressed); falling back to stream size");
      try( InputStream desc = am.open("interspec-assets.zip", AssetManager.ACCESS_RANDOM) )
      {
        long total = 0, n;
        while( (n = desc.skip(1024L * 1024 * 1024)) > 0 )
          total += n;
        bundledZipLen = total;
      }
    }
    Log.d("copyWtAssets", "interspec-assets.zip is " + bundledZipLen + " bytes");

    // Decide whether to extract: only if assets don't exist yet, or the bundled size differs
    // from the previously-recorded size.
    boolean needToExtract = !parent_dir.exists() || !sizeMarker.exists();
    if( !needToExtract )
    {
      try( Scanner reader = new Scanner(sizeMarker) )
      {
        final long prev = reader.nextLong();
        needToExtract = (prev != bundledZipLen);
        Log.d("copyWtAssets", "previous zip size=" + prev + ", current=" + bundledZipLen
                            + ", needToExtract=" + needToExtract);
      }
      catch( Exception e )
      {
        Log.w("copyWtAssets", "Could not read previous asset size, will re-extract", e);
        needToExtract = true;
      }
    }

    if( !needToExtract )
    {
      Log.d("copyWtAssets", "No need to extract; " + wtAssetsDir + " already current");
      return;
    }

    Log.d("copyWtAssets", "Extracting interspec-assets.zip into " + wtAssetsDir);

    // Recursive-delete the existing assets so files removed in a new version of the zip don't
    // linger.  The old code called File.delete() on the four top-level subdirs, which returns
    // false silently for non-empty dirs -- meaning the cleanup never happened in practice.
    if( parent_dir.exists() )
      deleteRecursive(parent_dir);
    parent_dir.mkdirs();

    final String parent_canonical = parent_dir.getCanonicalPath();
    final String parent_prefix = parent_canonical + File.separator;

    try( BufferedInputStream bis = new BufferedInputStream(am.open("interspec-assets.zip"));
         ZipInputStream zis = new ZipInputStream(bis) )
    {
      final byte[] buffer = new byte[64 * 1024];
      ZipEntry ze;
      while( (ze = zis.getNextEntry()) != null )
      {
        final File file = new File(parent_dir, ze.getName());
        // Zip-slip guard: ensure the resolved path stays under parent_canonical, comparing
        // against (parent_canonical + File.separator) so that a sibling like "wt-assets2" can't
        // satisfy startsWith("wt-assets") and escape.
        final String canonicalPath = file.getCanonicalPath();
        if( !canonicalPath.equals(parent_canonical) && !canonicalPath.startsWith(parent_prefix) )
          throw new SecurityException("Zip path traversal blocked: '" + canonicalPath
                                      + "' is outside '" + parent_canonical + "'");

        if( ze.isDirectory() )
        {
          file.mkdirs();
        }
        else
        {
          final File p = file.getParentFile();
          if( p != null ) p.mkdirs();
          try( FileOutputStream fos = new FileOutputStream(file) )
          {
            int count;
            while( (count = zis.read(buffer)) != -1 )
              fos.write(buffer, 0, count);
          }
        }
      }
    }

    Log.d("copyWtAssets", "Extraction complete");

    try( FileWriter wr = new FileWriter(sizeMarker) )
    {
      wr.write(Long.toString(bundledZipLen));
    }
    Log.d("copyWtAssets", "Wrote asset_size.txt");
  }

  /** Recursively deletes a file or directory and all its contents.  Returns true if everything
   *  was removed successfully.  Java 7 / API 21 compatible (no Files.walkFileTree). */
  private static boolean deleteRecursive( File f )
  {
    if( f.isDirectory() )
    {
      final File[] children = f.listFiles();
      if( children != null )
      {
        for( File child : children )
          deleteRecursive(child);
      }
    }
    return f.delete();
  }


  public static native void initNative();
  public static native int startServingInterSpec( String process_name, String userdatadir, String basedir, String xml_config_path );
  public static native int openFile( String sessionToken, String filepath);
  public static native boolean openAppUrl( String sessionToken, String url );
  public static native boolean killServer();
  public static native boolean setTempDir( String tmpdir );
  /** Java side previously declared this as (String require) while the JNI implementation took
   *  jboolean; calling it would have passed a jstring pointer where a jboolean was expected
   *  (undefined behavior).  Fixed signature: boolean in -> boolean out. */
  public static native boolean setRequireSessionToken( boolean require );
  public static native boolean addPrimarySessionToken( String token );
  public static native boolean addExternalSessionToken( String token );
  public static native int removeSessionToken( String token );
  public static native int setInitialFileToLoad( String token, String filepath );
  /** Routes a Wt-served download URL through InterSpecServer::download_to_native_save, which
   *  fetches the response in-process via Wt::Http::Client and pipes the bytes back through
   *  the native file-save handler registered in startServingInterSpec.  That handler in turn
   *  calls our CallbackFromNativeInterface, which launches the ACTION_CREATE_DOCUMENT save dialog. */
  public static native void nativeDownloadUrl( String url );
  static native void setFileSaveCallback(CallbackFromNativeInterface cb);
}