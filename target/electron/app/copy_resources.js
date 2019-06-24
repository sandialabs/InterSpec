var fs = require('fs.extra');
const path = require('path');

let dirsToCopy = [
  "C:/Users/wcjohns/Documents/InterSpec/build_electron/resources",
  "C:/Users/wcjohns/Documents/InterSpec/data",
  "C:/Users/wcjohns/Documents/InterSpec/InterSpec_resources",
  "C:/Users/wcjohns/Documents/InterSpec/example_spectra"
];

let executablesToCopy = ["C:/Users/wcjohns/Documents/InterSpec/build_electron/MinSizeRel/InterSpec.exe.exe"];

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

//Copy over SandiaDecay to the data directory - CMake probably already should have done this, but just to be sure.
var snldecay_src = "C:/Users/wcjohns/Documents/InterSpec/external_libs/SandiaDecay/sandia.decay.nocoinc.min.xml"
var snldecay_dest = "./data/sandia.decay.xml";

if( fs.lstatSync(snldecay_dest).isFile() ){
  fs.chmodSync(snldecay_dest, 0775);
  fs.unlinkSync(snldecay_dest);
}
fs.copy(snldecay_src, snldecay_dest, { replace: true }, function (err) {
  if (err)
    console.log("Failed to copy '" + snldecay_src + "' to '" + snldecay_dest + "'" + err);
  else
    console.log("Copied '" + snldecay_src + "' to '" + snldecay_dest + "'");
  
  fs.chmod(snldecay_dest, 0555, function(err){
    if(err) throw err;
    console.log( "Changed permissions of " + dest );
  });
});


//
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

