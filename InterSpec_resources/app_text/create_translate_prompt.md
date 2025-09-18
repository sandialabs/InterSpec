# Introduction

The XML files in the `InterSpec_resources/app_text` directory are used to provide the application text in the InterSpec application.  InterSpec is an application for performing gamma spectroscopy of radiation data.  The user audience of the application are highly technical and knowledgable scientists or Engineers.


# Information for translating
- English is the base, and "truth" information language.  All translations should be done from English, to the target laguage.
- The XML files given in English are named according to the ".cpp" file that they primarily provide strings for.  For example `InterSpec.xml` gives the strings for `src/InterSpec.cpp`, and `FluxTool.xml` gives strings for `src/FluxTool.cpp`.  The translated XML files all end with an underscore and the two-letter country code, and optionally a dash with specialization, and then the ".xml".  For example `InterSpec_fr.xml` contain the French strings, `InterSpec_de.xml` the german, `InterSpec_en-GB.xml` the British English strings, and so on.
- To add a language, you do not need to translate all elements in an XML file, or all XML files.  For any element, or file you do not translate, InterSpec will fall back to using the default (English) values for those strings.  
- Please make sure the translated file contains a "<message>" element for evenry element in the original English file.
- For strings that are less than or equal to about 5 words, you can assume that the string labels a GUI element (e.g., the label for an input, or the text in a button).  For these strings, translating in a compact and succinct manner is paramount, because if the string is too long, it may mess the layout of the application up.
- Translation of individual strings should take the context of the file and application into account, because this is how the user will be seeing the string - this may help to shorten some translations.
- Comments and aextra attributes in the XML should be kept in the translated XML.
- The English XML is frequently updated, so the translated XML files should be totally updated, when a new translation is requested.
- All elements of the English XML file should be in the translated XML file.  If for some reason a element cant be translated, or it is better to not translate it into the target language, the element should still be placed in the target XML file.
- The `id` attribute for each "<message>" XML element must remain identical, and not changed, as this attribute is what is used to match the string to where it goes in the application.
- White-space formatting (new lines, lines between elements), should be preserved when it can be, and when it makes sense.
- If you ask questions about the meaning of a string, please place a comment in the English XML providing the explanation or context, for next time.

# Instructions for translating

For each English XML file in `InterSpec_resources/app_text` please create a translated XML file, named according to the above rules, for the target language.  If anything is not clear, or there are questions about the meaning of a string - please see the InterSpec source code (your can search for the "id" attribute value in the source code), or ask questions!

Also, importantly, if the phrase being translated is under about 5 English words, please choose as short and succinct translation as you can (the string is probably labeling a GUI element, where space is often limited), making sure it still makes sense in context.

Also, please make sure to translate all English XML files to the target language, and translate each "<message>" entry in the original XML to the target XML.

# Language and Country specific instructions

None yet.