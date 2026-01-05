InterSpec.icon is the Icon Composer icon file.
InterSpec.icns and Assets.car are generated from the .icon file, using a command like:
```
xcrun actool InterSpec.icon --compile ./tmp --platform macosx --target-device mac --minimum-deployment-target 10.15 --app-icon InterSpec --include-all-app-icons --output-partial-info-plist /dev/null
```

The Assets.car file is for macOS Tahoe (26) and newer, while InterSpec.icns is for previous versions of macOS.