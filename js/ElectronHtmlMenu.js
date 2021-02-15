WT_DECLARE_WT_MEMBER
(TitleBarChangeMaximized, Wt::JavaScriptFunction, "TitleBarChangeMaximized",
function(maximized){
  if( maximized ){
    $('.resizer.top,.resizer.left').hide();
    $('.window-max-restore').removeClass('window-maximize').addClass('window-unmaximize');
  }else{
    $('.resizer.top,.resizer.left').show();
    $('.window-max-restore').removeClass('window-unmaximize').addClass('window-maximize');
  }
}
);


