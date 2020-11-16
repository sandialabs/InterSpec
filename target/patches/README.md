Process for patching files.

# To make a patch
```bash
cp someFile.txt someFile.txt.orig
```
Then edit someFile.txt, then

```bash
diff -Naur someFile.txt.orig  someFile.txt > someFile.txt.patch
```

# To apply a patch
```bash
patch -u someFile.txt -i someFile.txt.patch
```

# Command to apply patches for Wt 3.3.4
```bash
patch -u wt-3.3.4/CMakeLists.txt -i /path/to/InterSpec/target/patches/wt/3.3.4/CMakeLists.txt.patch
patch -u wt-3.3.4/src/Wt/Render/CssParser.C -i /path/to/InterSpec/target/patches/wt/3.3.4/CssParser.C.patch
patch -u wt-3.3.4/src/http/RequestParser.C -i /path/to/InterSpec/target/patches/wt/3.3.4/RequestParser.C.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr.C -i /path/to/InterSpec/target/patches/wt/3.3.4/ptr.C.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr.patch -i /path/to/InterSpec/target/patches/wt/3.3.4/ptr.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr_impl.h -i /path/to/InterSpec/target/patches/wt/3.3.4/ptr_impl.h.patch
patch -u wt-3.3.4/src/Wt/Dbo/Query.C -i /path/to/InterSpec/target/patches/wt/3.3.4/Query.C.patch
patch -u wt-3.3.4/src/Wt/WDllDefs.h -i /path/to/InterSpec/target/patches/wt/3.3.4/WDllDefs.h.patch
```