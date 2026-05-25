# ProGuard / R8 rules for InterSpec.
#
# minifyEnabled is currently false in app/build.gradle so these rules are not consulted at
# build time -- but the keep set has to be correct *before* anyone turns minification on,
# otherwise the first release build under R8 will silently strip JNI symbols and the JS
# bridge methods and the resulting APK will UnsatisfiedLinkError on startup.

# Keep the whole InterSpec Activity package by name so:
#  - JNI symbol resolution can find Java_gov_sandia_InterSpec_InterSpec_* methods, and
#  - the Activity class name embedded in AndroidManifest.xml still resolves at runtime.
-keep class gov.sandia.InterSpec.** { *; }
-keepnames class gov.sandia.InterSpec.**

# Native methods (the `native` keyword) must keep their original name/signature so the JVM
# can match them against the JNI exports in libInterSpecAppLib.so.
-keepclasseswithmembernames class * {
    native <methods>;
}

# @JavascriptInterface-annotated methods are called reflectively from JS running inside the
# WebView; if R8 renames or drops them the JS bridge stops working.
-keepclassmembers class * {
    @android.webkit.JavascriptInterface <methods>;
}

# Stack-trace line numbers in crash reports are much more useful than nothing; this keeps
# them without exposing the original source filename.
-keepattributes SourceFile,LineNumberTable
-renamesourcefileattribute SourceFile
