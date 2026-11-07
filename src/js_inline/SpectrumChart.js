/* Note: this is at the same time valid JavaScript and C++. */


WT_DECLARE_WT_MEMBER
(CanvasToPngDownload, Wt::JavaScriptFunction, "CanvasToPngDownload",
function(elid,filename)
{
  try{
    var canvas = document.getElementById(elid).querySelector('canvas');
    if( canvas.length === 0 )
      throw 'Could not find canvas';
      
    if( filename.length === 0 )
      filename = "spectrum.png"; //ToDo: add in date/time or something.
    var dt = canvas.toDataURL('image/png');
    dt = dt.replace(/^data:image\\/[^;]*/, 'data:application/octet-stream');
    dt = dt.replace(/^data:application\\/octet-stream/, 'data:application/octet-stream;headers=Content-Disposition%3A%20attachment%3B%20filename='+filename);
    
    var link = document.createElement("a");
    link.download = filename;
    link.href = dt;
    link.target="_blank";
    link.click();
  }catch(e){
    console.log( 'Error saving PNG from spectrum: ' + e );
  }
}
);

