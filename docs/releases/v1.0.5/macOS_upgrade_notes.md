# Upgrading from InterSpec 1.0.4 to 1.0.5 instructions on macOS
Starting in version 1.0.5 the InterSpec application bundle ID changed from "gov.sandia.InterSpec" to "gov.sandia.macOS.InterSpec".  Unfortunately this means macOS sees version 1.0.5 as a different application than 1.0.4, so none of your user preferences or work/spectra you've saved into `InterSpec`s internal database will be available after the upgrade.

If you don't care about retaining your previous preferences or work history, then no special action is needed when you upgrade.

If you would like to retain you previous preferences and work history, please complete the following steps
- Upgrade `InterSpec` by downloading the version 1.0.5 DMG file, and after opening, drag `InterSpec.app` to your `Applications` directory as normal when installing/upgrading an application
- Run the new version of `InterSpec`, and exit the application once it completes starting up.
- Open <b>Terminal.app</b> and copy-paste the following command into the terminal.
  ```bash
  cp ~/Library/Containers/sandia.InterSpec/Data/Library/Application\ Support/sandia.InterSpec/InterSpecUserData.db ~/Library/Containers/gov.sandia.macOS.InterSpec/Data/Library/Application\ Support/sandia.InterSpec/
  ```
- Restart `InterSpec` and all your old settings and work should now be restored within `InterSpec`


## Explanation
While preparing to have the app [notarized by Apple](https://developer.apple.com/documentation/security/notarizing_your_app_before_distribution), it was noted that the iOS and macOS versions of an app should not share the same bundle ID, as was previously the case for `InterSpec`.  
to prevent a conflict between macOS and iOS versions of the app.  
Starting in `InterSpec` version 1.0.5 the macOS InterSpec executable is , in addition to having the app sandbox and hardened runtime options enabled.
