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
 
package eu.webtoolkit.android.InterSpec;

import eu.webtoolkit.android.WtAndroid;
import android.app.Activity;
import android.content.res.AssetManager;
import android.os.Bundle;
import android.webkit.*;
import android.net.*;
import android.content.*;
import android.util.*;
import android.os.Handler;
import android.os.Message;
import android.view.View;
import android.os.Build;
import android.view.View;
import android.view.GestureDetector;
import android.view.MotionEvent;
import android.view.ViewGroup;
import android.provider.DocumentsProvider;
import android.os.ParcelFileDescriptor;
import java.io.FileDescriptor;
import java.nio.channels.FileChannel;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileOutputStream;
import android.provider.DocumentsContract;
import android.os.Environment;
import android.database.Cursor;
import android.provider.MediaStore;
import android.provider.OpenableColumns;
import java.util.Random;
import eu.webtoolkit.android.InterSpec.R;

	
public class InterSpec extends Activity
{
  private int httpPort = 0;
  private String interspecID = "";
  private View mDecorView;
  private ValueCallback<Uri> mUploadMessage;
  private final static int FILECHOOSER_RESULTCODE = 1;
 
 
  public final class WtWebChromeClient extends WebChromeClient 
  {
    @Override
    public boolean onJsAlert(WebView view, String url,  String message,  JsResult result) 
	{
      Log.d("WtAndroid::onJsAlert", message);
      return true;
    }

    //openFileChooser for Android < 3.0
	public void openFileChooser(ValueCallback<Uri> uploadMsg) 
	{
	  mUploadMessage = uploadMsg;
      Intent i = new Intent(Intent.ACTION_GET_CONTENT);
      i.addCategory(Intent.CATEGORY_OPENABLE);
      i.setType("*/*");
      InterSpec.this.startActivityForResult(
	    Intent.createChooser(i,"File Chooser"), FILECHOOSER_RESULTCODE);
    }
    
	//openFileChooser for Android 3.0+
    public void openFileChooser(ValueCallback<Uri> uploadMsg, String acceptType) 
	{
      openFileChooser(uploadMsg);
    }                   

    //openFileChooser for other Android versions
    public void openFileChooser(ValueCallback<Uri> uploadMsg, String acceptType, String capture) 
	{
      openFileChooser(uploadMsg);
    } 
	
    @Override
    public boolean onConsoleMessage (ConsoleMessage consoleMessage)
    {
      Log.d("WtWebChromeClient::onConsoleMessage", consoleMessage.message());
      return true;
    }
  }//public final class WtWebChromeClient extends WebChromeClient 


  public String getFilePath( final Uri uri ) 
  {
	//adapted from https://github.com/iPaulPro/aFileChooser/tree/master/aFileChooser
    final Context context = getApplicationContext();
	final boolean isKitKat = Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT;

    if( isKitKat && DocumentsContract.isDocumentUri(context, uri) ) 
	{
      // ExternalStorageProvider
      if( "com.android.externalstorage.documents".equals(uri.getAuthority()) ) 
	  {
        final String docId = DocumentsContract.getDocumentId(uri);
        final String[] split = docId.split(":");
        final String type = split[0];

        if( "primary".equalsIgnoreCase(type) )
		{
          return Environment.getExternalStorageDirectory() + "/" + split[1];
        }
        // TODO handle non-primary volumes
      }else if( "com.android.providers.downloads.documents".equals(uri.getAuthority()) ) 
	  {
        final String id = DocumentsContract.getDocumentId(uri);
        final Uri contentUri = ContentUris.withAppendedId(
            Uri.parse("content://downloads/public_downloads"), Long.valueOf(id));

        return getDataColumn(contentUri, null, null);
      } //else if media
    }else if ("content".equalsIgnoreCase(uri.getScheme())) 
	{
      return getDataColumn(uri, null, null);
    }else if ("file".equalsIgnoreCase(uri.getScheme())) 
	{
      return uri.getPath();
    }

    return null;
  }//public String getFilePath( final Uri uri ) 

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
	}
	
    return outputname;
  }//public void copyUriToTmpDir( Uri result ) 
  
  public void openFileInInterSpec( Uri result )
  {
    if( result == null )
      return;
	
    boolean shouldDeleteFile = false;
    String pathname = getFilePath( result );
    
    if( pathname == null )
    {
      shouldDeleteFile = true;
      String displayName = getDisplayFileName( result );	
      if( displayName == null )
        displayName = "unamedfile";
      pathname = copyUriToTmpDir( result, displayName );
    }//if( pathname == null )
		
    if( pathname != null )
    {
      //now we should tell InterSpec to open pathname.
      //  The only problem is we dont know if it should be foreground, sceond, or background
		  
	  Log.d("openFileInInterSpec", "Will send the following to InterSpec: " + pathname);
	  int staus = WtAndroid.openfileininterppec( pathname, 0, interspecID );  
	  
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
     if( null == mUploadMessage )
		 return;
	 
      Uri result = intent == null || resultCode != RESULT_OK 
		           ? null
                   : intent.getData();	  
	  
      openFileInInterSpec( result );
      
      mUploadMessage.onReceiveValue(result);
      mUploadMessage = null;
    }//if( pathname != null )
  }//protected void onActivityResult(...)
  

  @Override
  public void onCreate( Bundle savedInstanceState )
  {
	Log.d("onCreate", "Starting");
	  
	Random rn = new Random();
	interspecID = "androidsession" + Integer.toString(rn.nextInt(9999999));
		 
    if( httpPort == 0 )
    {
	  Log.d("onCreate", "Starting server");
      super.onCreate(savedInstanceState);
    
      this.requestWindowFeature( android.view.Window.FEATURE_NO_TITLE );
    
      WtAndroid.settmpdir( getCacheDir().getPath() );
    
      httpPort = WtAndroid.startWt(this);
    
      setContentView(R.layout.main);
        
      WebView webview = (WebView)findViewById(R.id.webview);
      webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID);
      
      webview.getSettings().setSupportMultipleWindows(false);
      webview.getSettings().setJavaScriptEnabled(true);
/*
      webview.getSettings().setUseWideViewPort(true);
      webview.getSettings().setAllowFileAccess(true);
*/
        
      webview.setWebChromeClient(new WtWebChromeClient());
	  
	  /*Enable remote debugging of webview, requires API 19 (KITKAT) */
	  if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT ) 
	  {
	    WebView.setWebContentsDebuggingEnabled(true);
	  }
    }//if( httpPort == 0 )

    Log.d("onCreate", "done starting server ish");
    
    Intent intent = getIntent();
    String action = intent.getAction();
    String type = intent.getType();
        
    if( Intent.ACTION_VIEW.equals(action) ) 
    {  
	  Log.d("onCreate", "ACTION_VIEW");
		
       Uri fileUri = (Uri) intent.getData();
       if( fileUri != null )
       {
          int entrynum = WtAndroid.addopenfiletodb( fileUri.getPath() );  
          WebView webview = (WebView)findViewById(R.id.webview);
          webview.loadUrl("http://localhost:" + httpPort + "/?externalid=" + interspecID + "&specfile=" + entrynum );
       }
    }else if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null ) 
    {
      // Need to handle multiple files here
	  Log.d("onCreate", "ACTION_SEND_MULTIPLE");
    }else 
    {
      // Handle other intents, such as being started from the home screen
	  Log.d("onCreate", "Another Action: " + action );
    }
    
	/* Need at least API 19 for the full screen emmersive view 
	 * I havent yet tested to make sure the bellow protections 
	 * for earlier API's actually works
     */
    if( Build.VERSION.SDK_INT >= Build.VERSION_CODES.KITKAT ) 
	{
	  mDecorView = getWindow().getDecorView();
      mDecorView.setOnSystemUiVisibilityChangeListener(
        new View.OnSystemUiVisibilityChangeListener() 
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
	   
    }else if( Intent.ACTION_SEND_MULTIPLE.equals(action) && type != null )
    {
        // Need to handle multiple files here
    } else
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
  
  
}