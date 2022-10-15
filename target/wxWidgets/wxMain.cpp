
#ifdef WX_PRECOMP
  #include "wx/wxprec.h"
#else
  #include "wx/wx.h"
#endif

#if !wxUSE_WEBVIEW_WEBKIT && !wxUSE_WEBVIEW_WEBKIT2 && !wxUSE_WEBVIEW_EDGE
#error "A wxWebView backend is required"
#endif






#include "InterSpecWxApp.h"


wxIMPLEMENT_APP(InterSpecWxApp);
