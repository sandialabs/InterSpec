#!/bin/bash
#
# Post-install script for the .deb and .rpm packages.
#
# This is a fork of electron-builder's own templates/linux/after-install.tpl, kept because
# `deb.afterInstall` / `rpm.afterInstall` *replace* that template rather than appending to it
# (app-builder-lib/out/targets/FpmTarget.js), so changing one line means owning the file.
#
# The only deliberate difference from upstream is the chrome-sandbox block below; everything else
# is upstream's and should be kept in sync with it.  The "Check electron-builder postinst template
# drift" step in .github/workflows/build_app.yml hashes the upstream template and fails the build
# when it changes, so improvements there (their AppArmor handling especially) do not silently pass
# us by.
#
# Note on templating: electron-builder substitutes ${executable} and ${sanitizedProductName} here,
# using the regex /\${([a-zA-Z]+)}/g and *throwing* on any macro it does not recognise.  So never
# write a brace-style shell variable whose name is only letters - a dollar-brace HOME or PATH
# would fail the whole build with "Macro HOME is not defined", even inside a comment like this one.
# Use the unbraced $HOME form, or names containing underscores, as upstream and the code below do.

if type update-alternatives >/dev/null 2>&1; then
    # Remove previous link if it doesn't use update-alternatives
    if [ -L '/usr/bin/${executable}' -a -e '/usr/bin/${executable}' -a "`readlink '/usr/bin/${executable}'`" != '/etc/alternatives/${executable}' ]; then
        rm -f '/usr/bin/${executable}'
    fi
    update-alternatives --install '/usr/bin/${executable}' '${executable}' '/opt/${sanitizedProductName}/${executable}' 100 || ln -sf '/opt/${sanitizedProductName}/${executable}' '/usr/bin/${executable}'
else
    ln -sf '/opt/${sanitizedProductName}/${executable}' '/usr/bin/${executable}'
fi

# --- BEGIN InterSpec deviation from upstream -------------------------------------------------
# Decide whether the SUID sandbox helper needs to be setuid-root.
#
# Upstream probes with a bare `unshare --user true`.  But this script runs as *root* during package
# configuration, while kernel.apparmor_restrict_unprivileged_userns (Ubuntu 23.10+) restricts only
# *unprivileged* namespace creation - root is exempt.  So upstream's probe always succeeds on those
# systems, concludes namespaces work, and leaves chrome-sandbox non-setuid for the unprivileged
# user who actually runs the app.  It is testing the wrong subject.
#
# That is masked wherever the AppArmor profile installed below loads, because the profile grants
# the app `userns` directly.  Where the profile cannot load - e.g. Ubuntu 23.10, whose AppArmor 3.x
# does not support our abi/4.0 profile, yet which does restrict unprivileged userns - the app ends
# up with neither mechanism and Chromium aborts during startup with
#   FATAL:sandbox/linux/suid/client/setuid_sandbox_host.cc: The SUID sandbox helper binary was
#   found, but is not configured correctly
# before any JavaScript is evaluated, so the app cannot recover from it itself.
#
# Probing as an unprivileged user tests the case that actually matters.  Erring toward setuid is
# safe: when the namespace sandbox does work, Chromium simply never uses the helper.
_interspec_unpriv_userns_works() {
  if command -v runuser >/dev/null 2>&1; then
    runuser -u nobody -- unshare --user true >/dev/null 2>&1
  elif command -v setpriv >/dev/null 2>&1; then
    setpriv --reuid=65534 --regid=65534 --clear-groups unshare --user true >/dev/null 2>&1
  else
    # Cannot drop privileges; fall back to upstream's probe rather than guessing.
    { [ -L /proc/self/ns/user ] && unshare --user true >/dev/null 2>&1; }
  fi
}

if ! _interspec_unpriv_userns_works; then
    # Unprivileged user namespaces are unavailable, so Chromium will need the SUID helper.
    chmod 4755 '/opt/${sanitizedProductName}/chrome-sandbox' || true
else
    chmod 0755 '/opt/${sanitizedProductName}/chrome-sandbox' || true
fi
# --- END InterSpec deviation from upstream ---------------------------------------------------

if hash update-mime-database 2>/dev/null; then
    update-mime-database /usr/share/mime || true
fi

if hash update-desktop-database 2>/dev/null; then
    update-desktop-database /usr/share/applications || true
fi

# Install apparmor profile. (Ubuntu 24+)
# First check if the version of AppArmor running on the device supports our profile.
# This is in order to keep backwards compatibility with Ubuntu 22.04 which does not support abi/4.0.
# In that case, we just skip installing the profile since the app runs fine without it on 22.04.
#
# Those apparmor_parser flags are akin to performing a dry run of loading a profile.
# https://wiki.debian.org/AppArmor/HowToUse#Dumping_profiles
#
# Unfortunately, at the moment AppArmor doesn't have a good story for backwards compatibility.
# https://askubuntu.com/questions/1517272/writing-a-backwards-compatible-apparmor-profile
if apparmor_status --enabled > /dev/null 2>&1; then
  APPARMOR_PROFILE_SOURCE='/opt/${sanitizedProductName}/resources/apparmor-profile'
  APPARMOR_PROFILE_TARGET='/etc/apparmor.d/${executable}'
  if apparmor_parser --skip-kernel-load --debug "$APPARMOR_PROFILE_SOURCE" > /dev/null 2>&1; then
    cp -f "$APPARMOR_PROFILE_SOURCE" "$APPARMOR_PROFILE_TARGET"

    # Updating the current AppArmor profile is not possible and probably not meaningful in a chroot'ed environment.
    # Use cases are for example environments where images for clients are maintained.
    # There, AppArmor might correctly be installed, but live updating makes no sense.
    if ! { [ -x '/usr/bin/ischroot' ] && /usr/bin/ischroot; } && hash apparmor_parser 2>/dev/null; then
      # Extra flags taken from dh_apparmor:
      # > By using '-W -T' we ensure that any abstraction updates are also pulled in.
      # https://wiki.debian.org/AppArmor/Contribute/FirstTimeProfileImport
      apparmor_parser --replace --write-cache --skip-read-cache "$APPARMOR_PROFILE_TARGET"
    fi
  else
    echo "Skipping the installation of the AppArmor profile as this version of AppArmor does not seem to support the bundled profile"
  fi
fi
