The text throughout InterSpec is loaded at runtime from the XML files in [InterSpec_resources/app_text](https://github.com/sandialabs/InterSpec/tree/master/InterSpec_resources/app_text).  These files are named according to the C++ file they contain the text for.  

That is `InterSpec.xml` and the `InterSpec_xx-YY.xml` files contain the strings that [src/InterSpec.cpp](https://github.com/sandialabs/InterSpec/blob/master/src/InterSpec.cpp) use.  The `InterSpec.xml` file contains the default (English) language strings, while `InterSpec_fr.xml` contain the French strings, `InterSpec_de.xml` the german, `InterSpec_en-GB.xml` the British English strings, and so on. 
 The country codes scheme used to augment the file name are given in https://datatracker.ietf.org/doc/html/rfc5646#appendix-A, with some more examples being "en", "en-US", "fr", "fr-FR", "es-ES", "ja" Japanese, "zh-Hant" Chinese written using the Traditional Chinese script, "ru" Russian, etc.


In these XML files, there are entries that look like:
```xml
  <message id="app-mi-file-about">About InterSpec</message>
  <message id="app-mi-file-open">Open File...</message>
  <message id="app-mi-file-loaded-spec">Loaded Spectra...</message>
  ...
```
The C++ code use the `id` attribute value to retrieve the `message` elements value.  
So to add a new language, you will copy the default (English) XML file(s) to new file names, with the "_xx-YY" of the 
new language appended to the XML file names, then you will translate the element for the `message` elements, while 
leaving the `id` attribute value unchanged.

For example, the French translation of the above three elements, would become:
```xml
  <message id="app-mi-file-about">À propos d'InterSpec</message>
  <message id="app-mi-file-open">Ouvrir un fichier...</message>
  <message id="app-mi-file-loaded-spec">Spectres chargés...</message>
```

***

In addition to the XML files in 
[InterSpec_resources/app_text](https://github.com/sandialabs/InterSpec/tree/master/InterSpec_resources/app_text),
you can also translate:
- [InterSpec_resources/static_text/help.json](https://github.com/sandialabs/InterSpec/tree/master/InterSpec_resources/static_text/help.json),
but for this file, you should only translate the `title` fields of this JSON, but use the same naming scheme for each language, as for the XML files.  This file defines the help content titles in the help dialog.
- The XML files in [InterSpec_resources/static_text/](https://github.com/sandialabs/InterSpec/tree/master/InterSpec_resources/static_text/), using the same naming scheme as for the other files.  These files are XHTML files that give the "help" content for each of the tools.


***

To add a language, you do not need to translate all elements in an XML file, or all XML files.  For any element, or file you do not translate, InterSpec will fall back to using the default (English) values for those strings.  

At a minimum, you must add a version of `InterSpec.xml` for the new language - preferably translating all elements, but this isn't required.  In the InterSpec app, the available languages (Help &rarr; Languages &rarr;) are populated according to which `InterSpec_xx-YY.xml` files are present.

But in addition `InterSpec.xml` to you really should, at a minimum, modify:
    - `D3SpectrumDisplayDiv.xml`, `CompactFileManager.xml`, `ReferencePhotopeakDisplay.xml`, `PeakModel.xml`, `PeakInfoDisplay.xml`, `EnergyCalTool.xml`, `IsotopeSearchByEnergy.xml`
This will translate the text seen in the application in its basic view, including all the tool-tabs that are open by default, so things will look consistent, until the user opens one of the specialized tools.


***

Some potential issues you may run into:
- The XML files must be valid XML.
- The list of HTML character entities that can be used is limited (and can be seen [here](https://github.com/emweb/wt/blob/b84925215d2b45879cf20c0cb340c4e7960d0c53/src/3rdparty/rapidxml/rapidxml_xhtml.hpp#L71)).  Instead you should directly use the UTF-8 characters.
- The XML files must be UTF-8 encoded.  On Windows, UTF-16 is commonly used by applications - please make sure UTF-8 encoding is used.
- Unfortunately, if the file is not valid XML, or a un-allowed character entity is used, the error-messaging is not great.  When the application is loaded, the file wont be used, but on Windows, no error messages are printed out.  On macOS or Linux, if you run the InterSpec executable from the command-line, an error will be printed out to stderr, giving the location of the error in the XML file.


***

A note on performance 20240521:
On an M1 mac, built in release mode, for local server, timing from the time the `InterSpec` class is created, to the ender of its first `render()` call (which I think should capture the added overhead of parsing the language XML, and localizing strings), I compared the pre-internationalization time, to the post-internationalization.
On the first load, of post-internationalized, the worst performance I got was 49 ms, however, was usually ~35 ms.  On pre-internationalized, typical was 33ms.  On second load (without restarting app), the timings were essentially the same between apps, with it being 10 ms for second load, and 6 ms for subsequent loads.
So all-in-all, internationalization does not cause a performance issue - assuming how I took the timings was valid.