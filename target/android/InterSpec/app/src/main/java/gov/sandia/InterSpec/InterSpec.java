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
import androidx.core.app.ActivityCompat;

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
import java.io.FileDescriptor;
import java.nio.channels.FileChannel;
import java.nio.ByteBuffer;
import java.util.Scanner;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileOutputStream;
import android.provider.DocumentsContract;
import android.os.Environment;
import android.database.Cursor;
import android.provider.OpenableColumns;
import android.content.res.AssetFileDescriptor;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import android.app.Fragment;
import java.lang.ref.WeakReference;

import eu.webtoolkit.android.WtAndroid;
//import gov.sandia.InterSpecMainActivity.InterSpecMainActivity;
//import gov.sandia.InterSpecMainActivity.InterSpecMainActivity.R;
//import android.R;

import gov.sandia.InterSpec.R;

//TODO: look at extending with AppCompatActivity instead of Activity.  WIll need to add entry into AndroidManifest.xml
//  To use a compatible theme, and also rpobably add a style to get rid of app bar.
public class InterSpec extends Activity
{
  private int httpPort = 0;
  private String interspecID = "";
  private View mDecorView;
  /** File upload callback for platform versions prior to Android 5.0 */
  private ValueCallback<Uri> mUploadMessage;
  /** File upload callback for Android 5.0+ */
  private ValueCallback<Uri[]> mUploadMessageMult;
  private final static int FILECHOOSER_RESULTCODE = 1;
 
 
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


  public String getDataColumn(Uri uri, String selection, String[] selectionArgs) 
  {
    //from https://github.com/iPaulPro/aFileChooser/tree/master/aFileChooser
    final Context context = getApplicationContext();
    Cursor cursor = null;
    final String column = "_data";
    final String[] projection = { column };

    try
	  {
      cursor = context.getContentResolver().query(uri, projection, selection, selectionArgs, null);
      if (cursor != null && cursor.moveToFirst()) 
	    {
        final int column_index = cursor.getColumnIndexOrThrow(column);
        return cursor.getString(column_index);
      }
    }finally 
    {
      if( cursor != null )
        cursor.close();
    }
	
    return null;
  }//public static String getDataColumn(Context context, Uri uri, String selection, String[] selectionArgs) 

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
	
    boolean shouldDeleteFile = true;

    String displayName = getDisplayFileName( result );
    Log.d("openFileInInterSpec", "displayName=" + displayName );

    if( displayName == null )
      displayName = "unamedfile";
    String pathname = copyUriToTmpDir( result, displayName );

		
    if( pathname != null )
    {
      //now we should tell InterSpec to open pathname.
      //  The only problem is we dont know if it should be foreground, sceond, or background
		  
	  Log.d("openFileInInterSpec", "Will send the following to InterSpec: " + pathname);
	  int staus = openfileininterppec( pathname, 0, interspecID );

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
    if( requestCode == FILECHOOSER_RESULTCODE )
    {
      if( null != mUploadMessage )
      {
        Log.d("onActivityRe", "null != mUploadMessage");

        Uri result = intent == null || resultCode != RESULT_OK ? null : intent.getData();

        openFileInInterSpec(result);

        mUploadMessage.onReceiveValue(result);
        mUploadMessage = null;
      }else if( mUploadMessageMult != null )
      {
        Log.d("onActivityRe", "mUploadMessageMult != null");

        Uri[] dataUris;
        if( intent != null )
        {
          try
          {
            dataUris = new Uri[]{Uri.parse(intent.getDataString())};
          } catch (Exception e)
          {
            Log.d( "onActivityRes", "Failed to parse dataUris" );
            dataUris = null;
          }

          Log.d( "onActivityRes", "dataUris.length=" + dataUris.length );
          for (int i = 0; i < dataUris.length; i++)
          {
            openFileInInterSpec( dataUris[i] );
          }

          mUploadMessageMult.onReceiveValue(dataUris);
        }else
        {
          mUploadMessageMult.onReceiveValue( null );
        }//if( intent != null )


       mUploadMessageMult = null;
     }else
     {
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

    Random rn = new Random();
	  interspecID = "androidsession" + Integer.toString(rn.nextInt(9999999));
		 
    if( httpPort == 0 )
    {
	    Log.d("onCreate", "Starting server");
      super.onCreate(savedInstanceState);
    
      this.requestWindowFeature( android.view.Window.FEATURE_NO_TITLE );

      settmpdir( getCacheDir().getPath() );
    
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
      if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT) {
          webview.setLayerType(View.LAYER_TYPE_HARDWARE, null);
      } else {
          webview.setLayerType(View.LAYER_TYPE_SOFTWARE, null);
      }

      /*
      webview.getSettings().setUseWideViewPort(true);
      webview.getSettings().setAllowFileAccess(true);
      */

      webview.setWebChromeClient(new WtWebChromeClient());


      webview.setDownloadListener(new DownloadListener() {
        @Override
        public void onDownloadStart(final String url, final String userAgent, String contentDisposition, String mimetype, long contentLength)
        {
          Log.d("onDownloadStart", "Entering onDownloadStart");
          //checking runtime permissions
          if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.M) {
            if (checkSelfPermission(android.Manifest.permission.WRITE_EXTERNAL_STORAGE)
                    == PackageManager.PERMISSION_GRANTED) {

              downloadDialog(url,userAgent,contentDisposition,mimetype);
            } else {

              //requesting permissions
              ActivityCompat.requestPermissions(InterSpec.this, new String[]{Manifest.permission.WRITE_EXTERNAL_STORAGE}, 1);
            }
          }
          else {
            //Code for devices below API 23 or Marshmallow
            downloadDialog(url,userAgent,contentDisposition,mimetype);
          }

        }
      });

      if( Intent.ACTION_VIEW.equals(action) )
      {
        Log.d("onCreate", "ACTION_VIEW: Will load a file via URL");
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
              int entrynum = addopenfiletodb(pathname);
              Log.d("onCreate", "ACTION_VIEW, after copying to temp file '" + pathname
                      + "' will load from from dbentry=" + entrynum );
              webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID + "&specfile=" + entrynum);

              //We should delete the file at some point...
            }else
            {
              Log.d("onCreate", "ACTION_VIEW, Failed to copy resource to a temporary file." );
              webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID );
            }
          }else
          {
            // handle as file uri
            Log.d("onCreate", "ACTION_VIEW, fileUri is NOT a SCHEME_CONTENT, but a path=" + fileUri.getPath() );
            int entrynum = addopenfiletodb( fileUri.getPath() );
            Log.d("onCreate", "ACTION_VIEW, will load from from dbentry=" + entrynum );
            webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID + "&specfile=" + entrynum );
          }
        }else
        {
          Log.d("onCreate", "ACTION_VIEW, but not fileUri, will call loadUrl though");
          webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID );
        }
      }else
      {
        if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null )
        {
          Log.d("onCreate", "ACTION_SEND_MULTIPLE"); // Need to handle multiple files here
        }else
        {
          Log.d("onCreate", "Another Action: " + action );
        }

        Log.d("onCreate", "Loading default initial page");
        webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID);
      }

	    /*Enable remote debugging of webview, requires API 19 (KITKAT) */
	    if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT )
	    {
	      WebView.setWebContentsDebuggingEnabled(true);
	    }
    }//if( httpPort == 0 )

    Log.d("onCreate", "done starting server ish");

	/* Need at least API 19 for the full screen emmersive view */
    if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT ) 
	  {
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
	}//if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT ) 
  }//public void onCreate(Bundle savedInstanceState)


  public void downloadDialog(final String url,final String userAgent,String contentDisposition,String mimetype)
  {
    //filename using url.
    final String filename = URLUtil.guessFileName(url,contentDisposition,mimetype);
    //creates alertdialog
    AlertDialog.Builder builder=new AlertDialog.Builder(InterSpec.this);
    //alertdialog title
    builder.setTitle("Download");
    //alertdialog message
    builder.setMessage("Do you want to save " +filename);

    builder.setPositiveButton("Yes", new DialogInterface.OnClickListener()
    {
      @Override
      public void onClick(DialogInterface dialog, int which)
      {
        //DownloadManager.Request created with url.
        DownloadManager.Request request = new DownloadManager.Request(Uri.parse(url));
        //cookie
        String cookie=CookieManager.getInstance().getCookie(url);
        //Add cookie and User-Agent to request
        request.addRequestHeader("Cookie",cookie);
        request.addRequestHeader("User-Agent",userAgent);
        //file scanned by MediaScannar
        request.allowScanningByMediaScanner();
        //Download is visible and its progress, after completion too.
        request.setNotificationVisibility(DownloadManager.Request.VISIBILITY_VISIBLE_NOTIFY_COMPLETED);
        //DownloadManager created
        DownloadManager downloadManager=(DownloadManager)getSystemService(DOWNLOAD_SERVICE);
        //saves file in Download folder
        request.setDestinationInExternalPublicDir(Environment.DIRECTORY_DOWNLOADS, filename);
        //download enqued
        downloadManager.enqueue(request);
      }
    });
    builder.setNegativeButton("Cancel", new DialogInterface.OnClickListener() {
      @Override
      public void onClick(DialogInterface dialog, int which)
      {
        //cancel the dialog if Cancel clicks
        dialog.cancel();
      }

    });
    //alertdialog shows
    builder.create().show();

  }


  @Override
  protected void onNewIntent( Intent intent )
  {
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
    }
    else if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null )
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
	
	if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT ) 
    {
	  mHideHandler.removeMessages(0);
      if( hasFocus ) 
      {
        mHideHandler.sendEmptyMessageDelayed(0, 300);
      }
    }//if( KITKAT or greater )
  }//public void onWindowFocusChanged(boolean hasFocus) 
  
  
  private void hideSystemUI() 
  {
	if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT ) 
	{
      mDecorView.setSystemUiVisibility(
                  View.SYSTEM_UI_FLAG_LAYOUT_HIDE_NAVIGATION
                  | View.SYSTEM_UI_FLAG_LAYOUT_FULLSCREEN
                  | View.SYSTEM_UI_FLAG_HIDE_NAVIGATION
                  | View.SYSTEM_UI_FLAG_FULLSCREEN
                  | View.SYSTEM_UI_FLAG_LOW_PROFILE
                  | View.SYSTEM_UI_FLAG_IMMERSIVE_STICKY);
	//View.SYSTEM_UI_FLAG_IMMERSIVE View.SYSTEM_UI_FLAG_IMMERSIVE_STICKY
    }
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
    Log.d("WtAndroid::onCreate", "in startWt ...");

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

    List<String> args = new ArrayList<String>();
    args.add("app");

    args.add("--docroot");
    args.add(wtAssetsDir);
    args.add("--approot");
    args.add(wtAssetsDir);

    if( userDataDir.isEmpty() )
    {
      args.add("--userdatadir");
      args.add(userDataDir);
    }

    args.add("--http-address");
    args.add("127.0.0.1");
    args.add("--http-port");
    args.add("0");

    args.add("-c");
    args.add(wtAssetsDir + "/data/config/wt_config_android.xml");

    File tmpDir = new File(activity.getFilesDir().getAbsolutePath() + "/tmp");
    tmpDir.mkdir();

    //ToDo: delete tmpDir contents if it exists
    args.add("-DWT_TMP_DIR=" + tmpDir.getAbsolutePath());

    //Quite down all the log outputs of ever GET or REQUEST
    args.add("--accesslog=-");

    //Electron version of app also specifies "--userdatadir", "--basedir", and "--externalid"

    Log.d("WtAndroid::onCreate", "Starting wt application ...");
    String[] argv = new String[args.size()];
    args.toArray(argv);


    int httpPort = WtAndroid.startwt(argv);
    Log.d("WtAndroid::onCreate", "Started wt application on http-port " + httpPort);

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


  public static native int addopenfiletodb( String path );
  public static native int openfileininterppec( String path, int type, String sessionid );
  public static native void settmpdir( String tmpPath );
}