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

import android.Manifest;
import android.app.Activity;
import android.app.DownloadManager;
import android.content.pm.PackageManager;
import android.content.DialogInterface;
import android.content.res.AssetManager;
import android.os.Bundle;
import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;

import android.provider.DocumentsContract;
import android.provider.MediaStore;
import android.webkit.*;
import android.net.*;
import android.content.*;
import android.util.*;
import android.os.Handler;
import android.os.Message;
import android.view.View;
import android.os.Build;
import android.view.GestureDetector;
import android.view.MotionEvent;
import android.os.ParcelFileDescriptor;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileDescriptor;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.List;
import java.util.Scanner;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileOutputStream;
import android.os.Environment;
import android.database.Cursor;
import android.provider.OpenableColumns;
import android.widget.Toast;

import java.io.InputStream;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.concurrent.Executors;

//import gov.sandia.InterSpecMainActivity.InterSpecMainActivity;
//import gov.sandia.InterSpecMainActivity.InterSpecMainActivity.R;
//import android.R;

import gov.sandia.InterSpec.R;

interface CallbackInterface {
  public void callback();
}

public class InterSpec extends AppCompatActivity
{
  private int httpPort = 0;
  private String interspecID = "";
  private View mDecorView;
  /** File upload callback for platform versions prior to Android 5.0 */
  private ValueCallback<Uri> mUploadMessage;
  /** File upload callback for Android 5.0+ */
  private ValueCallback<Uri[]> mUploadMessageMult;
  private final static int FILECHOOSER_RESULTCODE = 1;

  private final static int SAVE_AS_CODE = 9821341;

  // We're keeping \c mUrlToDownload around for the moment, in hopes we can undo the Android specific
  //  hack for downloading files
  private String mUrlToDownload;

  // Using hacked saving to temporary file in Android, instead of via network download of file.
  //   These next two variables keep track of the temporary  file made by the c++ we will need to
  //   copy over to where the user selects.
  private String mTempFileToSave;
  private String mTempFileDisplayName;


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

    @Override
    public boolean onConsoleMessage (ConsoleMessage consoleMessage)
    {
      Log.d("onConsoleMessage", consoleMessage.message());
      return true;
    }
  }//public final class WtWebChromeClient extends WebChromeClient 


  public static void saveFileFromTempLocation( String location, String displayName )
  {
    // This function can probably get deleted See TRY_TO_CALL_JAVA_FROM_CPP in android.cpp for notes
    Log.d("fromCpp", "Being called from C++ location=" + location + ", diaplyName=" + displayName );
  }

  public String getDisplayFileName( Uri result ) 
  {
	  String name = null;
    Cursor cursor = InterSpec.this.getContentResolver().query(result, null, null, null, null, null);
    try 
    {
	  if( cursor != null && cursor.moveToFirst() ) 
        name = cursor.getString( cursor.getColumnIndex(OpenableColumns.DISPLAY_NAME) );
    }finally
    {
	  if( cursor != null )
        cursor.close();
    }
	
	  return name;
  }//public String getDisplayFileName( Uri result ) 

  public String copyUriToTmpDir( Uri result, final String displayName ) 
  {
	  String outputname = null;
  	try
  	{
      ParcelFileDescriptor parcelFD = InterSpec.this.getContentResolver().openFileDescriptor(result, "r");
      FileDescriptor fd = parcelFD.getFileDescriptor();
      FileChannel inchannel = new FileInputStream(fd).getChannel();
	    outputname = InterSpec.this.getFilesDir().getAbsolutePath()
		                        + java.io.File.separator + "tmp" 
							    + java.io.File.separator + displayName;
	
      File outputfile = new File( outputname );
      outputfile.createNewFile();
      FileChannel outchannel = new FileOutputStream(outputfile).getChannel();
	
      try 
      {
        inchannel.transferTo(0, inchannel.size(), outchannel);
      }finally
      {
        if( inchannel != null )
         inchannel.close();
        if( outchannel != null )
          outchannel.close();
      }
    }catch( java.io.IOException e )
	  {
      Log.d("copyUriToTmpDir", "IOException copying file" );
	    return null;
	  }catch( java.lang.NullPointerException e )  //possibly thrown by the parcelFD.getFileDescriptor() line
    {
      Log.d("copyUriToTmpDir", "NullPointerException" );
      return null;
    }
	
    return outputname;
  }//public void copyUriToTmpDir( Uri result ) 
  
  public void openFileInInterSpec( Uri result )
  {
    Log.d("openFileInInterSpec", "Got URI" );

    if( result == null )
      return;

    if( result.getScheme().contentEquals("interspec") )
    {
      Log.d("openFileInInterSpec", "Got scheme of interspec, url: " + result.toString()  );
      openAppUrl( interspecID, result.toString() );
      return;
    }


    boolean shouldDeleteFile = true;

    String displayName = getDisplayFileName( result );
    Log.d("openFileInInterSpec", "displayName=" + displayName );

    if( displayName == null )
      displayName = "unamedfile";
    String pathname = copyUriToTmpDir( result, displayName );

		
    if( pathname != null )
    {
      //now we should tell InterSpec to open pathname.
      //  The only problem is we dont know if it should be foreground, secondary, or background
		  
	  Log.d("openFileInInterSpec", "Will send the following to InterSpec: " + pathname);
	  // TODO: as long as there are no errors, we can clean up the file immediately, but if there is an error we should do something...
      int status = openFile( interspecID, pathname );

      if( status <= 0 )
      {
        Log.d("onActivityResult", "Failed to open file with status: " + status );
      }

      Log.d("onActivityResult", "Open file status: " + pathname);

      if( shouldDeleteFile )
      {
        File f = new File( pathname );	
        boolean deleted = f.delete();
        if( !deleted )
        {
          Log.d("onActivityResult", "Couldnt cleaned up: " + pathname);
        }
	  }//if( shouldDeleteFile )
	}//if( pathname != null )
  }//public void openFileInInterSpec( Uri result )

  @Override
  protected void onActivityResult( int requestCode, int resultCode, Intent intent )
  {
    super.onActivityResult(requestCode, resultCode, intent);

    if( (requestCode == SAVE_AS_CODE) ){
      // Using hacked saving to temporary file in Android, instead of via network download of file.

      if( resultCode != Activity.RESULT_OK ) {
        if( !mTempFileToSave.isEmpty() ){
          File file = new File( mTempFileToSave );
          if( file.delete() )
            Log.d("DD run", "Deleted temporary file: " + mTempFileToSave );
          else
            Log.d("DD run", "Failed to delete temporary file: " + mTempFileToSave );
        }// if( we have a temporary file we need to cleanup )
        return;
      }//if( user decided to not save file )

      final Uri uri = intent.getData();

       // Downloading the link via the DownloadManager, or the below manual connection always fails;
       //  we get zero bytes.
       //  I'm guessing the failure is due to Wt's session management; I would think this would be
       //  either using cookies, or the session-id; however, the session-id is already in the
       //  download URL, and we are setting the cookie - so I really dont know why it isnt working.
       //  I basically gave up and am hacking in having the c++ copy files to temporary files, and
       //  then copying those to the URI the user selects.

      /*
       // This is code I would *think* should work to copy the data from a download URL; leaving in
       //  for future investigations.
       // Using hacked saving to temporary file in Android, instead of via network download of file.
       WebView webview = (WebView)findViewById(R.id.webview);
       final String appUrl = webview.getUrl();
       final String url = mUrlToDownload;
       String cookie = CookieManager.getInstance().getCookie(url);
       String userAgent = webview.getSettings().getUserAgentString();

      ExecutorService myExecutor = Executors.newCachedThreadPool();
      myExecutor.execute(() -> {
        Log.d("DD run", "Will download: " + url + " to " + uri.getPath() );

        Context context = getApplicationContext();
        try
        {
          OutputStream output = context.getContentResolver().openOutputStream(uri);

          HttpURLConnection con = (HttpURLConnection) new URL(url).openConnection();
          con.setRequestMethod("GET");
          con.addRequestProperty("Cookie", cookie );
          con.addRequestProperty("User-Agent", userAgent );
          BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));

          //All of this is crap - I cant figure out how to download a local url - its amazing
          //  Next step is to really start thinking terribly and differently.

          // Diff buffer types...I to keep things happy; I dont know.
          byte dataBuffer[] = new byte[8*1024];
          char charBuffer[] = new char[8*1024];

          int bytesRead, totalBytes = 0;

          while ((bytesRead = in.read(charBuffer, 0, 8*1024)) != -1) {
            for( int i = 0; i < bytesRead; ++i )
              dataBuffer[i] = (byte)charBuffer[i];
            output.write(dataBuffer, 0, bytesRead);
            totalBytes += bytesRead;
          }

          output.flush();
          output.close();
          Log.d("DD run", "Done writing from URL: " + url + " - " + totalBytes + " total bytes." );

          //String filename = uri.getLastPathSegment();
          //Log.d("DD run", "Will make: " + filename + " visible to user." );

          //Make it so user can actually see the file in the Downloads folder
          //DownloadManager downloadManager=(DownloadManager)getSystemService(DOWNLOAD_SERVICE);
          //downloadManager.addCompletedDownload( uri.getPath(), "Spectrum File", true, "application/octet-stream", filename, totalBytes, true);

        } catch (IOException e) {
          Toast.makeText(context, "Error writing output file", Toast.LENGTH_SHORT).show();
        }
      });
*/

      // Using hacked saving to temporary file in Android, instead of via network download of file.
      Log.d("DD run", "Got location: " + mTempFileToSave + " and displayName:" + mTempFileDisplayName );

      Context context = getApplicationContext();

      if( mTempFileToSave.isEmpty() )
      {
        Toast.makeText(context, "No output file to copy to destination - sorry!", Toast.LENGTH_SHORT).show();
        return;
      }

      try
      {
        // Copy file to URI, then delete temp file
        OutputStream output = context.getContentResolver().openOutputStream(uri);
        FileInputStream is = new FileInputStream(mTempFileToSave);

        byte[] buffer = new byte[1024];
        int length, totalBytes = 0;
        while ((length = is.read(buffer)) > 0) {
          output.write(buffer, 0, length);
          totalBytes += length;
        }

        output.flush();
        output.close();
        is.close();

        Log.d("DD run", "Done writing: " + mTempFileToSave + " - " + totalBytes + " total bytes - to dest" );

        // Right now the file wont show up in the Downloads app; we would have to call
        //  downloadManager.addCompletedDownload(...) below, but we dont actually have the filename
        //  or path, and it appears nigh impossible to get this from the URI (i.e., the URI may
        //  not even point to a file on the filesystem), so I guess we just wont notify the
        //  download manager (and the DownloadManager only takes in "file" URIs).
        //DownloadManager downloadManager=(DownloadManager)getSystemService(DOWNLOAD_SERVICE);
        //downloadManager.addCompletedDownload( uri.getPath(), "Spectrum File", true, "application/octet-stream", mTempFileDisplayName, totalBytes, true);

        Toast.makeText(context, "Done writing file to disk.", Toast.LENGTH_SHORT).show();
      } catch (IOException e) {
        Toast.makeText(context, "Error writing output file", Toast.LENGTH_SHORT).show();
      }

      if( !mTempFileToSave.isEmpty() )
      {
        File file = new File( mTempFileToSave );
        if (file.delete()) {
          Log.d("DD run", "Deleted temporary file: " + mTempFileToSave );
        }
        else {
          Log.d("DD run", "Failed to delete temporary file: " + mTempFileToSave );
        }
      }//if( !mTempFileToSave.isEmpty() )

      mTempFileDisplayName = "";
      mTempFileToSave = "";

      return;
    }//if( this is result of save as dialog )

    if (requestCode == FILECHOOSER_RESULTCODE) {
      if (null != mUploadMessage) {
        Log.d("onActivityRe", "null != mUploadMessage");

        Uri result = intent == null || resultCode != RESULT_OK ? null : intent.getData();

        openFileInInterSpec(result);

        mUploadMessage.onReceiveValue(result);
        mUploadMessage = null;
      } else if (mUploadMessageMult != null) {
        Log.d("onActivityRe", "mUploadMessageMult != null");

        Uri[] dataUris;
        if (intent != null) {
          try {
            dataUris = new Uri[]{Uri.parse(intent.getDataString())};
          } catch (Exception e) {
            Log.d("onActivityRes", "Failed to parse dataUris");
            dataUris = null;
          }

          Log.d("onActivityRes", "dataUris.length=" + dataUris.length);
          for (int i = 0; i < dataUris.length; i++) {
            openFileInInterSpec(dataUris[i]);
          }

          mUploadMessageMult.onReceiveValue(dataUris);
        } else {
          mUploadMessageMult.onReceiveValue(null);
        }//if( intent != null )


        mUploadMessageMult = null;
      } else {
        Log.d("onActivityRe", "Neither upload messages were non - null - unexpected!");
      }
    }//if( pathname != null )
  }//protected void onActivityResult(...)

  @Override
  public void onCreate( Bundle savedInstanceState )
  {
    Log.d("onCreate", "Starting");

    Intent intent = getIntent();
    String action = intent.getAction();
    String type = intent.getType();
    Uri data = intent.getData();

    if( data != null )
      Log.d("onCreate", "data: " + data.toString() );
    else
      Log.d("onCreate", "data: null" );

    Random rn = new Random();
    interspecID = "session" + Integer.toString(rn.nextInt(9999999));

    addPrimarySessionToken( interspecID );

    if( httpPort == 0 )
    {
      initNative();

      setFileSaveCallback( new CallbackInterface() {
        public void callback() {
          // TODO: in the future we could maybe have this callback (or one that takes a cople string arguments)
          //       launch the Save As dialog, and then trigger the copying of files
          Log.d("callback", "Being called-back when file is to be saved");
        }
      });

      setTempDir( getCacheDir().getPath() );

      Log.d("onCreate", "Starting server");
      super.onCreate(savedInstanceState);
    
      this.requestWindowFeature( android.view.Window.FEATURE_NO_TITLE );
      httpPort = startWt(this);
    
      setContentView(R.layout.main);
        
      WebView webview = (WebView)findViewById(R.id.webview);
      WebSettings settings = webview.getSettings();
      settings.setCacheMode(WebSettings.LOAD_NO_CACHE);
      settings.setSupportMultipleWindows(false);
      settings.setJavaScriptEnabled(true);
      settings.setLoadWithOverviewMode(true);
      settings.setUseWideViewPort(true);
      settings.setSupportZoom(true);
      settings.setGeolocationEnabled(false);
      settings.setBuiltInZoomControls(false);
      settings.setLayoutAlgorithm(WebSettings.LayoutAlgorithm.SINGLE_COLUMN);
      settings.setDomStorageEnabled(true);
      webview.setScrollBarStyle(WebView.SCROLLBARS_OUTSIDE_OVERLAY);
      webview.setScrollbarFadingEnabled(true);
      webview.setLayerType(View.LAYER_TYPE_HARDWARE, null);

      /*
      webview.getSettings().setUseWideViewPort(true);
      webview.getSettings().setAllowFileAccess(true);
      */

      webview.setWebChromeClient( new WtWebChromeClient() );
/*
      // Using hacked saving to temporary file in Android, instead of via network download of file.
      DownloadListener ourDownloadListner = new DownloadListener() {
        @Override
        public void onDownloadStart(final String url, final String userAgent, String contentDisposition, String mimetype, long contentLength)
        {
          Log.d("onDownloadStart", "Entering onDownloadStart: url=" + url + ", Length: " + contentLength );
          ... have user select location ... and get the data ...
      };//new DownloadListener(){...}

      webview.setDownloadListener( ourDownloadListner );
*/


      // TODO: I think see https://stackoverflow.com/questions/3926629/downloadlistener-not-working to get download listner working
      WebViewClient ourWebClient = new WebViewClient(){

        // you tell the webclient you want to catch when a url is about to load
        @Override
        public boolean shouldOverrideUrlLoading(WebView  view, String  url){
          Log.d("shouldOverrideUrlLoad", "shouldOverrideUrlLoading: " + url );

          // Using hacked saving to temporary file in Android, instead of via network download of file.
          //  Once we use the callback (set via setFileSaveCallback) to trigger saving of file (after
          //  some more testing), then this function of overriding URL
          String[] spoolAndDisplay = mostRecentSaveLocation();
          if( spoolAndDisplay.length == 0 )
          {
            Toast.makeText( getApplicationContext(), "Error writing output file - perhaps Android workaround for saving files isnt implemented in the C++ for this case.", Toast.LENGTH_SHORT).show();
            return true;
          }

          String location = spoolAndDisplay[0];
          String displayName = spoolAndDisplay[1];
          mTempFileToSave = location;
          mTempFileDisplayName = displayName;
          Log.d("shouldOverrideUrlLoad", "Got location: " + location + " and displayName:" + displayName );

          Intent intent = new Intent(Intent.ACTION_CREATE_DOCUMENT);
          intent.addCategory(Intent.CATEGORY_OPENABLE);
          intent.setType("*/*"); //not needed, but maybe usefull
          intent.putExtra(Intent.EXTRA_TITLE, displayName); //not needed, but maybe usefull
          mUrlToDownload = url;
          startActivityForResult(intent, SAVE_AS_CODE);

/*
          // Using hacked saving to temporary file in Android, instead of via network download of file.
          //view.loadUrl(url);
          //view.setDownloadListener(ourDownloadListner);
          ... let user select where to save ... and download file
 */
          return true;
        }
      };
      webview.setWebViewClient( ourWebClient );



      if( Intent.ACTION_VIEW.equals(action) )
      {
        Log.d("onCreate", "ACTION_VIEW: Will set initial file to load at start");
        Uri fileUri = (Uri) intent.getData();

        if( fileUri != null )
        {
          String scheme = fileUri.getScheme();
          if( ContentResolver.SCHEME_CONTENT.equals(scheme) )
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
        Log.d("onCreate", "ACTION_SEND_MULTIPLE"); // Need to handle multiple files here
      }else
      {
        Log.d("onCreate", "Another Action: " + action );
      }

      Log.d("onCreate", "Will load 'http://127.0.0.1:" + httpPort + "/?apptoken=" + interspecID + "'");
      webview.loadUrl("http://127.0.0.1:" + httpPort + "/?apptoken=" + interspecID );

      WebView.setWebContentsDebuggingEnabled(true);
    }//if( httpPort == 0 )

    Log.d("onCreate", "done starting server ish");

    mDecorView = getWindow().getDecorView();
    mDecorView.setOnSystemUiVisibilityChangeListener( new View.OnSystemUiVisibilityChangeListener()
    {
        @Override
        public void onSystemUiVisibilityChange(int flags)
        {
        boolean visible = (flags & View.SYSTEM_UI_FLAG_HIDE_NAVIGATION) == 0;
//			WebView webview = (WebView)findViewById(R.id.webview);
        if( visible )
        {
          //make it so we are no longer in emmersive mode, so this way
          //  Android will resize the WebView so that the active form
          //  (assuming keyboard is visible) will be visible.  There
          //  is an issue that when you hide the keyboard, the system
          //  navigation UI will still be visible until you tap somewhere
          //  else (whach calls hideSystemUI()).
          //  Non-optimal, but I cant figure out how to detect if the
          //  keyboard is visible.
          mDecorView.setSystemUiVisibility(0);
        }
      }//public void onSystemUiVisibilityChange(int flags)
    });
	 
    final View contentView = findViewById( R.id.webview );
    contentView.setClickable( true );
    final GestureDetector clickDetector = new GestureDetector( this,
            new GestureDetector.SimpleOnGestureListener()
    {
          @Override
          public boolean onSingleTapUp(MotionEvent e) 
          {
            boolean visible = (mDecorView.getSystemUiVisibility()
	                          & View.SYSTEM_UI_FLAG_HIDE_NAVIGATION) == 0;
            if (visible) 
            {
              hideSystemUI();
            } 
            return false;
          }
    });
	   
    contentView.setOnTouchListener(
        new View.OnTouchListener() 
        {
          @Override
          public boolean onTouch(View view, MotionEvent motionEvent) 
          {
            return clickDetector.onTouchEvent(motionEvent);
          }
    });
  }//public void onCreate(Bundle savedInstanceState)

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

    if( Intent.ACTION_VIEW.equals(action) )
    {
      Uri fileUri = (Uri) intent.getData();

      openFileInInterSpec( fileUri );
    } else if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null )
    {
        // Need to handle multiple files here
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
    super.onDestroy();
  }
  
  
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
  
  private final Handler mHideHandler = new Handler() 
  {
    @Override
    public void handleMessage(Message msg) 
	{
      hideSystemUI();
    }
  };
  

	  
  
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


  private static void copyWtAssets(String wtAssetsDir, AssetManager am) throws IOException
  {
    //We need to know if the assets have been updated since we last ran things to see if we need to
    // re-extract things (takes ~6 seconds to extract).  So as a hack for doing this, we will see
    // if interspect-assets.zip has changed size since last time, and only extract it if it has.
    // I know this isnt quite right, but close enough for now.

    boolean needToExtract = !(new File(wtAssetsDir).exists());

    if( !needToExtract )
    {
      Log.d("copyWtAssets", "initially wtAssetsDir exists");
      needToExtract = !(new File(wtAssetsDir + "/asset_size.txt").exists());
      Log.d("copyWtAssets", "then needToExtract=" + needToExtract );
    }

    if( !needToExtract )
    {
      Log.d("copyWtAssets", "In !needToExtract" );

      try
      {
        Log.d("copyWtAssets", "Will check size of interspect-assets.zip");

        InputStream desc = am.open("interspec-assets.zip", AssetManager.ACCESS_RANDOM );

        //Only takes ~1/10 second to get the filesize this way
        long ziplen = 0, nread = 0;
        while( (nread = desc.skip(1024*1024*1024)) != 0 ) {
          ziplen += nread;
        }
        Log.d("copyWtAssets", "ziplen = " + ziplen + " bytes");

        //Takes ~1/10 second to read in the previous file size.
        Scanner reader = new Scanner(new File(wtAssetsDir + "/asset_size.txt"));
        long i = reader.nextLong();
        Log.d("copyWtAssets", "previous ziplen = " + i + " bytes" );
        needToExtract = (i != ziplen);
      }catch( Exception e )
      {
        Log.d("copyWtAssets", "Caught exception trying to compare sizes" );
        needToExtract = true;
      }
    }


    if( needToExtract )
    {
      Log.d("copyWtAssets", "We do need to extract" );

      if( new File(wtAssetsDir).exists() ) {
        try{
          new File(wtAssetsDir + "/data").delete();
          new File(wtAssetsDir + "/InterSpec_resources").delete();
          new File(wtAssetsDir + "/resources").delete();
          new File(wtAssetsDir + "/example_spectra").delete();
        }catch( Exception e ){
          //Only happens when security manager deinfined - i guess.
          Log.d("copyWtAssets", "Caught exception trying delete previous resource directory" );
        }
      }else {
        new File(wtAssetsDir).mkdir();
      }

      long ziplen = 0;
      {
        InputStream desc = am.open("interspec-assets.zip", AssetManager.ACCESS_RANDOM );
        long nread = 0;

        while( (nread = desc.skip(1024*1024*1024)) != 0 ) {
          ziplen += nread;
        }

        Log.d("copyWtAssets", "interspec-assets.zip declared length is " + ziplen + " bytes" );
      }

      BufferedInputStream bis = new BufferedInputStream(am.open("interspec-assets.zip"));

      ZipInputStream zis = new ZipInputStream(bis);
      try
      {
        byte[] buffer = new byte[1024];
        int count;

        ZipEntry ze;

        while ((ze = zis.getNextEntry()) != null)
        {
          String file = wtAssetsDir + ze.getName();
          if (ze.isDirectory())
          {
            new File(file).mkdirs();
          } else
          {
            FileOutputStream fos = new FileOutputStream(file);
            while ((count = zis.read(buffer)) != -1)
            {
              fos.write(buffer, 0, count);
            }
            fos.close();
          }
        }
      } finally
      {
        zis.close();
      }

      Log.d("copyWtAssets", "Finished extracting interspec-assets.zip");

      FileWriter wr = new FileWriter(wtAssetsDir + "/asset_size.txt");
      wr.write(ziplen + "");
      wr.close();

      Log.d("copyWtAssets", "Wrote asset_size.txt");
    } else
    {
      Log.d("copyWtAssets", "No need to extract interspec-assets.zip, " + wtAssetsDir + " already exists" );
    }
  }


  public static native void initNative();
  public static native int startServingInterSpec( String process_name, String userdatadir, String basedir, String xml_config_path );
  public static native int openFile( String sessionToken, String filepath);
  public static native boolean openAppUrl( String sessionToken, String url );
  public static native boolean killServer();
  public static native boolean setTempDir( String tmpdir );
  public static native boolean setRequireSessionToken( String require );
  public static native boolean addPrimarySessionToken( String token );
  public static native boolean addExternalSessionToken( String token );
  public static native int removeSessionToken( String token );
  public static native int setInitialFileToLoad( String token, String filepath );
  public static native String[] mostRecentSaveLocation();
  static native void setFileSaveCallback(CallbackInterface cb);
}