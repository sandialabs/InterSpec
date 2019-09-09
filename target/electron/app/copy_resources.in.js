var fs = require('fs.extra');
const path = require('path');

let dirsToCopy = [
  "${PROJECT_BINARY_DIR}/resources",
  "${PROJECT_SOURCE_DIR}/data",
  "${PROJECT_SOURCE_DIR}/InterSpec_resources",
  "${PROJECT_SOURCE_DIR}/example_spectra"
];

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
var snldecay_src = "${PROJECT_SOURCE_DIR}/external_libs/SandiaDecay/sandia.decay.nocoinc.min.xml"
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

