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
 
package eu.webtoolkit.android;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import android.app.Activity;
import android.content.res.AssetManager;
import android.os.Bundle;
import android.util.Log;
import android.view.WindowManager;
import android.webkit.JsResult;
import android.webkit.WebChromeClient;
import android.webkit.WebView;


public class WtAndroid extends Activity
{	
  public static int startWt(Activity activity)
    {	
	    Log.d("WtAndroid::onCreate", "Extracting wt-assets.zip ...");
        
	String wtAssetsDir 
		  =activity.getFilesDir().getAbsolutePath()  + "/wt-assets/";
	try {
	    		copyWtAssets(wtAssetsDir, activity.getAssets());
	} catch (IOException e) {
			e.printStackTrace();
	}
	
	Log.d("WtAndroid::onCreate", "Finished extracting wt-assets.zip");
		
	List<String> args = new ArrayList<String>();
	args.add("app");
        
	args.add("--docroot");
	args.add(wtAssetsDir);
	args.add("--approot");
	args.add(wtAssetsDir);
        
	args.add("--http-address");
	args.add("127.0.0.1");
	args.add("--http-port");
	args.add("0");

	args.add("-c");
	args.add(wtAssetsDir + "/data/config/wt_config_android.xml");
 
	File tmpDir = new File(activity.getFilesDir().getAbsolutePath() + "/tmp");
	tmpDir.mkdir();
	args.add("-DWT_TMP_DIR=" + tmpDir.getAbsolutePath());
        
	Log.d("WtAndroid::onCreate", "Starting wt application ...");
	String[] argv = new String[args.size()];
        args.toArray(argv);
    	int httpPort = startwt(argv);
    	Log.d("WtAndroid::onCreate", "Started wt application on http-port " + httpPort);
	
	return httpPort;
    }
    
  private static void copyWtAssets(String wtAssetsDir, AssetManager am) throws IOException 
  {
    if( !new File(wtAssetsDir).exists() )
    {
      new File(wtAssetsDir).mkdir();
    	
      BufferedInputStream bis = new BufferedInputStream(am.open("wt-assets.zip"));
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
    }
  }


    private static native int startwt(String[] argv);
    public static native int addopenfiletodb( String path );
    public static native int openfileininterppec( String path, int type, String sessionid );
	  public static native void settmpdir( String tmpPath );
    
    static {
        System.loadLibrary("InterSpec");
    }
}
