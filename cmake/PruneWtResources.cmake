# PruneWtResources.cmake
#
# Copy Wt's resources/ directory into the build tree while removing the files
# InterSpec does not use, so packaged apps (iOS, Android, macOS, wxWidgets,
# Electron, web) bundle ~5 MB less of a ~5.9 MB tree.
#
# The removable set (WT_RESOURCES_PRUNE_LIST) is NOT hand-written: it is
# generated and verified by target/tools/wt_resource_audit/wt_resource_audit.py, which traces how
# Wt loads each resource and whether the current InterSpec source ever triggers
# that path.  See target/tools/wt_resource_audit/wt_resources_manifest.json for the per-file evidence, and
# re-run the script after a Wt upgrade (it also has a --check drift gate).
#
# Style mirrors cmake/DeployJsAndCss.cmake.

# Pull in the generated, evidence-backed removal list (sets WT_RESOURCES_PRUNE_LIST).
set( _wt_prune_list_file "${CMAKE_CURRENT_LIST_DIR}/wt_resources_prune_list.cmake" )
if( EXISTS "${_wt_prune_list_file}" )
  include( "${_wt_prune_list_file}" )
else()
  message( WARNING "PruneWtResources: ${_wt_prune_list_file} missing -- run "
                   "target/tools/wt_resource_audit/wt_resource_audit.py.  Will copy the full Wt "
                   "resources without pruning." )
  set( WT_RESOURCES_PRUNE_LIST )
endif()


# prune_wt_resources( <src_dir> <dst_dir> )
#
# Copies the *contents* of <src_dir> into <dst_dir> (basename-independent, so the
# destination is exactly <dst_dir> whatever the source dir is named), then removes
# every WT_RESOURCES_PRUNE_LIST entry from <dst_dir>.  <src_dir> is never touched,
# so passing a read-only / shared Wt install prefix is safe.
function( prune_wt_resources src_dir dst_dir )
  message( "Pruning unused Wt resources: ${src_dir} -> ${dst_dir}" )

  # Fresh copy each configure so source changes propagate and previously-removed
  # files never linger.
  file( REMOVE_RECURSE "${dst_dir}" )
  file( MAKE_DIRECTORY "${dst_dir}" )
  file( GLOB _items LIST_DIRECTORIES true "${src_dir}/*" )
  if( _items )
    file( COPY ${_items} DESTINATION "${dst_dir}" )
  endif()

  foreach( _rel ${WT_RESOURCES_PRUNE_LIST} )
    if( EXISTS "${dst_dir}/${_rel}" )
      file( REMOVE_RECURSE "${dst_dir}/${_rel}" )
    else()
      # Staleness alarm: a listed path absent from this Wt version means the
      # manifest is out of date (e.g. after a Wt upgrade that renamed it).
      message( WARNING "PruneWtResources: prune entry '${_rel}' not present in "
                       "${src_dir} -- stale list?  Re-run "
                       "target/tools/wt_resource_audit/wt_resource_audit.py." )
    endif()
  endforeach()
endfunction()
