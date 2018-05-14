var fs = require('fs.extra');
const path = require('path');

let dirsToCopy = [
  "${PROJECT_BINARY_DIR}/resources",
  "${PROJECT_SOURCE_DIR}/data",
  "${PROJECT_SOURCE_DIR}/InterSpec_resources",
  "${PROJECT_SOURCE_DIR}/example_spectra"
];

let executablesToCopy = ["${INTERSPEC_EXE_LOCATION}"];

for( dir of dirsToCopy ){
  
  console.log( "Will copy '" + dir + "' to '" + path.join('.',path.basename(dir)) + "'");
  
  fs.copyRecursive(dir, path.join('.',path.basename(dir)), function (err) {
    if (err) {
      console.log("Failed to copy '" + dir + "' to '.'" + err);
      //throw err;
    }
    console.log("Copied '" + dir + "' to '.'");
  });
}

for( file of executablesToCopy ){
  var dest = path.join('.',path.basename(file));
  
  //on windows MSVS will create InterSpec.exe.exe (due to my lame CMakeLists.txt) - make this sane
  dest = dest.replace( ".exe.exe", ".exe" );
  
  console.log( "Will copy '" + file + "' to '" + dest + "'");
  
  try{
    if( fs.lstatSync(dest).isFile() ){
      fs.chmodSync(dest, 0775);
      fs.unlinkSync(dest);
      console.log( "Removed old '" + dest + "'");
    }
  }catch( e ){
    console.log( "Failed to remove old '" + dest + "'");
  }
  

  fs.copy(file, dest, { replace: true }, function (err) {
    if (err) {
      console.log("Failed to copy '" + file + "' to '.'" + err);
      //throw err;
    }
    console.log("Copied '" + file + "' to '" + dest + "'");
    
    fs.chmod(dest, 0555, function(err){
      if(err) throw err;
      console.log( "Changed permissions of " + dest );
    });
  });
}

