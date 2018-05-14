#!/usr/bin/ruby

# InterSpec: an application to analyze spectral gamma radiation data.
# 
# Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
# For questions contact William Johnson via email at wcjohns@sandia.gov, or
# alternative emails of interspec@sandia.gov, or srb@sandia.gov.
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


###################################################################
#
# SCRIPT TO BUILD InterSpec LIBRARIES
#
###################################################################

require 'fileutils'

$REDIRECT = true

#------------------------------------------------------------------
#
# Class:       ElapsedTime
# Description: a class to report elapsed time
#
#------------------------------------------------------------------
class ElapsedTime

    #------------------------------------------------------------------
    #
    # Method:      initialize
    # Description: Constructor
    #
    #------------------------------------------------------------------
    def initialize
        @start_time = Time.new
        @elapsed = 0
    end

    #------------------------------------------------------------------
    #
    # Method:      format_rounded
    # Description: Round a number and format it using the given units
    # Param:       val - the number to format
    # Param:       units - either "second" or "minute"
    # Return:      a string ("1 second", "3 minutes", etc)
    #
    #------------------------------------------------------------------
    def format_rounded(val, units)
        rounded = val.round
        if rounded == 1
            return "1 #{units}"
        else
            return "#{rounded} #{units}s"
        end
    end

    #------------------------------------------------------------------
    #
    # Method:      format
    # Description: Format the elapsed time since this object was created
    # Return:      a string of the form "<n> seconds" or "<n> minutes and <n> seconds"
    #
    #------------------------------------------------------------------
    def format
        @elapsed = Time.now - @start_time
        seconds = format_rounded(@elapsed % 60, "second")
        minutes = format_rounded((@elapsed / 60).floor, "minute")
        if (@elapsed < 60)
            return seconds
        else
            return "#{minutes} and #{seconds}"
        end
    end
end

#------------------------------------------------------------------
#
# Class:       PreFlight
# Description: a class to check build conditions
#
#------------------------------------------------------------------
class PreFlight

    #------------------------------------------------------------------
    #
    # Method:      run_command
    # Description: run a shell command
    # Param:       cmd - the command to run
    # Return:      the output of the command as one big string
    #
    #------------------------------------------------------------------
    def run_command(cmd)
        output = `#{cmd}`
        return output
    end

    #------------------------------------------------------------------
    #
    # Method:      abort
    # Description: Print an error message and exit this script
    # Param:       msg - the text to print
    # Return:      none
    #
    #------------------------------------------------------------------
    def abort(msg)
        puts "    Aborting: #{msg}"
        exit 1
    end

    #------------------------------------------------------------------
    #
    # Method:      banner
    # Description: Display a message inside some banner text
    # Param:       msg - the message to display
    # Return:      none
    #
    #------------------------------------------------------------------
    def banner(msg)
        puts ''
        puts "----------------------------------------------------------------"
        puts "- #{msg}"
        puts "----------------------------------------------------------------"
    end


    #------------------------------------------------------------------
    #
    # Method:      get_xcode_version
    # Description: Return the XCode version number
    # Return:      a string like "4.5.2"
    #
    #------------------------------------------------------------------
    def get_xcode_version
        output = run_command("xcodebuild -version")
        # look for XCode n.n.n
        #puts output
        output =~ /XCode\s+(\d+\.\d+\.\d+)/i
        return $1
    end

    #------------------------------------------------------------------
    #
    # Method:      check_xcode
    # Description: Ensure that the required version of XCode is installed
    # Param:       required_version - the minimum XCode version number
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_xcode(required_version)
        banner("Checking XCode version")
        ver = get_xcode_version
        puts "    Your Xcode version is: #{ver}"
        abort("Need later version of Xcode") if ver < required_version
        puts "    Xcode version OK"
    end


    #------------------------------------------------------------------
    #
    # Method:      check_sdk
    # Description: Ensure that a particular SDK is installed
    # Param:       required - the sdk name (e.g. "iphoneos6.1")
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_sdk(required)
        result = false
        banner("Checking SDK versions")
        output = `xcodebuild -version -sdk`
        lines = output.split(/\n/)
        print "    You have the following SDKs: "
        lines.each do |line|
            if line =~ /^iPhoneOS.+\((.+)\)/ then
                print "#{$1} "
                result = $1 == required
            end
            if line =~ /^iPhoneSimulator.+\((.+)\)/ then
                print " #{$1}"
            end
        end
        print "\n";
        print("Did not find SDK version #{required}") unless result
        puts "    SDK version OK (#{required})"
    end

    #------------------------------------------------------------------
    #
    # Method:      check_directory
    # Description: Make sure a directory exists
    # Param:       dir - path to directory
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_directory(dir)
        abort("Required directory '#{dir}' does not exist") unless File.directory?(dir)
    end

    #------------------------------------------------------------------
    #
    # Method:      check_file
    # Description: Make sure a file exists
    # Param:       dir - path to the file
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_file(file)
        abort("Required file '#{file}' does not exist") unless File.file?(file)
    end


    #------------------------------------------------------------------
    #
    # Method:      check_directories
    # Description: Ensure that all build directories exist
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_directories
        banner("Checking required directories and files")
        check_directory('./target/ios/boost')
        check_file('./target/ios/boost/boost.sh')
        check_directory('./target/ios/wt/wt-3.3.4/src')
        check_file('./target/ios/wt/build-framework.sh')
        check_directory('/etc/wt')
        puts("    All required directories and files exist")
    end

end


#------------------------------------------------------------------
#
# Class:       Builder
# Description: a class perform InterSpec build operations
#
#------------------------------------------------------------------
class Builder


    #------------------------------------------------------------------
    #
    # Method:      initialize
    # Description: Constructor
    #
    #------------------------------------------------------------------
    def initialize
        @boost_dir = './target/ios/boost'
        @boost_libs = [
            "libboost_chrono.a",
            "libboost_date_time.a",
            "libboost_exception.a",
            "libboost_filesystem.a",
            "libboost_iostreams.a",
            "libboost_program_options.a",
            "libboost_random.a",
            "libboost_regex.a",
            "libboost_serialization.a",
            "libboost_signals.a",
            "libboost_system.a",
            "libboost_thread.a",
            "libboost_timer.a",
            "libboost_atomic.a"
        ]
        @wt_dir = './target/ios/wt/wt-3.3.4'

    end

    #------------------------------------------------------------------
    #
    # Method:      run_cmd
    # Description: Echo a command and then run it
    # Param:       cmd - the command to run
    # Return:      none
    #
    #------------------------------------------------------------------
    def run_cmd(cmd)
        puts("    Running: #{cmd}")
        system(cmd)
    end


    #------------------------------------------------------------------
    #
    # Method:      chdir
    # Description: Change directories
    # Param:       path - target directory
    #
    #------------------------------------------------------------------
    def chdir(path)
        puts("    chdir to #{path}")
        Dir.chdir(path)
    end


    def mkdir(path)
        if File.directory?(path)
            puts("    Directory #{path} already exists")
        else
            puts("    Creating directory #{path}")
            FileUtils.mkdir_p(path)
        end
    end

    def copyfile(from, to)
        puts("    Copy #{from} to #{to}")
        FileUtils.cp_r(from, to)
    end

    def rename(from, to)
        puts("    Rename #{from} to #{to}")
        FileUtils.mv(from, to)
    end

    #------------------------------------------------------------------
    #
    # Method:      abort
    # Description: Print an error message and exit this script
    # Param:       msg - the text to print
    # Return:      none
    #
    #------------------------------------------------------------------
    def abort(msg)
        puts "    Aborting: #{msg}"
        exit 1
    end



    #------------------------------------------------------------------
    #
    # Method:      check_directory
    # Description: Make sure a directory exists
    # Param:       dir - path to directory
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_directory(dir)
        abort("Required directory '#{dir}' does not exist") unless File.directory?(dir)
    end

    #------------------------------------------------------------------
    #
    # Method:      check_file
    # Description: Make sure a file exists
    # Param:       dir - path to the file
    # Return:      none (aborts on failure)
    #
    #------------------------------------------------------------------
    def check_file(file)
        abort("File '#{file}' does not exist") unless File.file?(file)
    end

    #------------------------------------------------------------------
    #
    # Method:      banner
    # Description: Display a message inside some banner text
    # Param:       msg - the message to display
    # Param:       pause_seconds - number of seconds to pause after banner
    # Return:      none
    #
    #------------------------------------------------------------------
    def banner(msg, pause_seconds = 0)
        puts ''
        puts "----------------------------------------------------------------"
        if msg.kind_of?(Array)
            msg.each do |m|
                puts "- #{m}"
            end
        else
            puts "- #{msg}"
        end
        puts "----------------------------------------------------------------"
        sleep(pause_seconds) unless pause_seconds == 0
    end

    #------------------------------------------------------------------
    #
    # Method:      delete_directory
    # Description: Delete a directory if it exists
    # Param:       path - the path to the directory
    # Return:      none (aborts if we are unable to delete the directory)
    #
    #------------------------------------------------------------------
    def delete_directory(path)
        if File.directory?(path) then
            puts("    Deleting directory: #{path}")
            FileUtils.rm_rf(path)
            abort("#{path} was not deleted") if File.directory?(path)
        end
    end

    #------------------------------------------------------------------
    #
    # Method:      delete_file
    # Description: Delete a file if it exists
    # Param:       path - the path to the file
    # Return:      none
    #
    #------------------------------------------------------------------
    def delete_file(path)
        if File.file?(path) then
            puts("    Deleting file: #{path}")
            File.delete(path)
        end
    end

    #------------------------------------------------------------------
    #
    # Method:      clean_boost
    # Description: Delete Boost build output directories
    # Return:      none (aborts if we can't delete one of the directories)
    #
    #------------------------------------------------------------------
    def clean_boost
        banner("Cleaning boost build products")
        delete_directory(@boost_dir + "/build")
        delete_directory(@boost_dir + "/framework")
        delete_directory(@boost_dir + "/prefix")
        delete_directory(@boost_dir + "/src")
    end

    #------------------------------------------------------------------
    #
    # Method:      clean_wt
    # Description: Delete WT build output directories
    # Return:      none (aborts if we can't delete one of the directories)
    #
    #------------------------------------------------------------------
    def clean_wt
        banner("Cleaning WT build products")
        delete_directory(@wt_dir + "/build")
        delete_directory(@wt_dir + "/build-armv7")
    end

    #------------------------------------------------------------------
    #
    # Method:      build_boost
    # Description: Build the boost libraries
    # Return:      none
    #
    #------------------------------------------------------------------
    def build_boost
        banner("Building boost in #{@boost_dir}")
        root_dir = Dir.pwd
        elapsed = ElapsedTime.new
        File.delete('./boost_build.log') if File.file?('./boost_build.log')
        if ($REDIRECT)
            puts("    Build output will be in boost_build.log")
            Dir.chdir(@boost_dir)
            run_cmd("IPHONEOS_DEPLOYMENT_TARGET=\"10.2\" ./boost.sh > ../../../boost_build.log  2>&1 ")
        else
            Dir.chdir(@boost_dir)
            run_cmd('./boost.sh')
        end
        etime = elapsed.format
        puts("    Boost build completed in #{etime}")
        Dir.chdir(root_dir)

#        rename("./target/ios/boost/ios/prefix", "./target/ios/prefix")
# copy( "target/ios/build/libs/boost/lib/allarch/*", "./target/ios/prefix/lib/" )

        @boost_libs.each do |f| check_file("./target/ios/prefix/lib/#{f}") end
    end

    #------------------------------------------------------------------
    #
    # Method:      build_wt
    # Description: Build the WT libraries
    # Return:      none
    #
    #------------------------------------------------------------------
    def build_wt
        build_dir = "target/ios/wt/build"
        root_dir = Dir.pwd
        banner("Building WT in #{build_dir}")
        abort("    The WT directory, #{@wt_dir}, does not exist") unless File.directory?("#{@wt_dir}")
        abort("    /etc/wt does not exist") unless File.directory?("/etc/wt")
        copyfile("cmake/WtFindBoost.txt", "#{@wt_dir}/cmake")
        copyfile("cmake/WtFindBoost-vintage.txt", "#{@wt_dir}/cmake")
        Dir.mkdir(build_dir) unless File.directory?(build_dir)
        logfile = 'wt_build.log'
        File.delete(logfile) if File.file?(logfile)

#As of 3.3.4 at least, WCanvasPaintDevice.C no longer needs patching
#        patchfile = "target/ios/wt/wt_3.3.1_WCanvasPaintDevice.c.patch";
#        patchtarget = "#{@wt_dir}/src/Wt/WCanvasPaintDevice.C";
#        check_file(patchfile)
#        check_file(patchtarget)
#        run_cmd("patch '--backup --forward' #{patchtarget} < #{patchfile}")
#        puts("Patched file ")

        puts("    Build output will be in #{logfile}")
        elapsed = ElapsedTime.new
        chdir(build_dir)
        script = "../build-framework.sh"
        check_file(script)
        if ($REDIRECT)
            run_cmd("IPHONEOS_DEPLOYMENT_TARGET=\"10.2\" bash #{script} > ../../wt_build.log 2>&1")
        else
            run_cmd("IPHONEOS_DEPLOYMENT_TARGET=\"10.2\" bash #{script}")
        end
        etime = elapsed.format
        puts("    WT build completed in #{etime}")
        # cp lib/*.a ~/install/ios/prefix/lib/
        chdir(root_dir)
    end

    #------------------------------------------------------------------
    #
    # Method:      copy_wt_files
    # Description: Copy WT library and include files to the boost build directory
    # Return:      none
    #
    #------------------------------------------------------------------
    def copy_wt_files
        build_dir = "target/ios/wt/build"
        root_dir = Dir.pwd
        prefix_dir = "#{root_dir}/target/ios/boost/ios/prefix"
        banner("Copying WT files into #{prefix_dir}")
        Dir.chdir(build_dir)
        run_cmd("cp stage/lib/*.a #{prefix_dir}/lib/")
        run_cmd("cp -r stage/include/Wt #{prefix_dir}/include/")
        run_cmd("cp -r stage/share #{prefix_dir}/")
        run_cmd("cp -r stage/cmake #{prefix_dir}/")
        chdir(root_dir)
    end

    #------------------------------------------------------------------
    #
    # Method:      rebuild_boost
    # Description: Clean the boost build outputs and rebuild it
    # Return:      none
    #
    #------------------------------------------------------------------
    def rebuild_boost
        clean_boost
        build_boost
    end


    #------------------------------------------------------------------
    #
    # Method:      rebuild_wt
    # Description: Clean the boost WT outputs and rebuild it
    # Return:      none
    #
    #------------------------------------------------------------------
    def rebuild_wt
        clean_wt
        build_wt
    end



    #------------------------------------------------------------------
    #
    # Method:      build_interspec
    # Description: Build the InterSpec library files for IOS
    # Return:      none
    #
    #------------------------------------------------------------------
    def build_interspec
        replaced_cmakelists = false
        banner(["Build the InterSpec library files for IOS"], 1.0)
        elapsed = ElapsedTime.new
        root_dir = Dir.pwd
        run_cmd('rm -rf build_ios')
        Dir.mkdir('build_ios')
        Dir.chdir('build_ios')
        run_cmd('which g++')
        run_cmd('which gcc')
        run_cmd("IPHONEOS_DEPLOYMENT_TARGET=\"10.2\" cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cmake/iOSToolchain.cmake -DIOS=ON ..")
        if ($REDIRECT)
            print "Build output is being redirected to ../interspec_build.log\n"
            run_cmd("IPHONEOS_DEPLOYMENT_TARGET=\"10.2\" make -j8 > ../interspec_build.log 2>&1")
        else
            run_cmd("IPHONEOS_DEPLOYMENT_TARGET=\"10.2\" make -j8")
        end
        chdir(root_dir)
        if (replaced_cmakelists) then
            copyfile('CMakeLists.save', 'CMakeLists.txt')
        end
        etime = elapsed.format
        puts("    InterSpec build completed in #{etime}")
    end


    def create_xcode_directories
        banner("Building directories for Xcode")
        mkdir './XcodeFiles/example_spectra'
        mkdir './XcodeFiles/data'
        mkdir './XcodeFiles/InterSpec_resources'
        mkdir './XcodeFiles/WtsRsrcs'
        mkdir './target/ios/include'
        run_cmd('cp -r ./data/* ./XcodeFiles/data')
        run_cmd('cp -r ./example_spectra/* ./XcodeFiles/example_spectra')
        run_cmd('cp -r ./Interspec_resources/* ./XcodeFiles/InterSpec_resources')
        run_cmd('cp -r ./target/ios/prefix/share/Wt/resources/* ./XcodeFiles/WtsRsrcs')
        run_cmd('rm  ./XcodeFiles/InterSpec_resources/static_text/videos/*ogv' )
        run_cmd('rm  ./XcodeFiles/data/wt_config.xml' )
#        run_cmd('rm  ./XcodeFiles/data/em_xs_data/*.xs' )
#        run_cmd('rm  ./XcodeFiles/WtsRsrcs/themes/polished' )
#        run_cmd('rm  ./XcodeFiles/WtsRsrcs/themes/bootstrap' )

#        run_cmd('cp -r ./target/ios/boost/ios/prefix/include/boost ./target/ios/include')
#        run_cmd('cp -r ./wt-3.3.4/build/stage/include/Wt ./target/ios/include')
#        run_cmd('cp -r ./InterSpec ./target/ios/include')
    end

    #------------------------------------------------------------------
    #
    # Method:      clean_interspec
    # Description: Delete the build_ios directory where InterSpec is built
    # Return:      none
    #
    #------------------------------------------------------------------
    def clean_interspec
        banner("Cleaning InterSpec build directory (build_ios)")
        delete_directory('build_ios')
    end



    #------------------------------------------------------------------
    #
    # Method:      build_all
    # Description: Build all the necessary IOS libraries
    # Return:      none
    #
    #------------------------------------------------------------------
    def build_all
        elapsed = ElapsedTime.new
        clean_all
        build_boost
        build_wt
        build_interspec
        create_xcode_directories
        etime = elapsed.format
        puts("BuildAll completed in #{etime}")
    end

    #------------------------------------------------------------------
    #
    # Method:      clean_all
    # Description: Clean the build outputs for Boost, WT, and InterSpec
    # Return:      none
    #
    #------------------------------------------------------------------
    def clean_all
        clean_boost
        clean_wt
        clean_interspec
    end


end


def contains(str, substring)
    return (str.index(substring) != nil)
end

#------------------------------------------------------------------
#
# Method:      banner
# Description: Display a message inside some banner text
# Param:       msg - the message to display
# Param:       pause_seconds - number of seconds to pause after banner
# Return:      none
#
#------------------------------------------------------------------
def banner(msg, pause_seconds = 0)
    puts ''
    puts "----------------------------------------------------------------"
    if msg.kind_of?(Array)
        msg.each do |m|
            puts "- #{m}"
        end
    else
        puts "- #{msg}"
    end
    puts "----------------------------------------------------------------"
    sleep(pause_seconds) unless pause_seconds == 0
end


def show_usage
    banner(["Use: build-ios [configtest | boost | wt | interspec | XcodeDirs | clean | all]",
        "Options:",
        "    boost     build the boost libraries",
        "    wt        build wt",
        "    interspec build InterSpec",
        "    XcodeDirs create directories to hold stuff needed by Xcode",
        "    all       build everything",
        "    clean     remove build products",
        " ",
        "You can specify more than one option on the command line,",
        " eg: ./iosbuild clean boost wt",
        " ",
        "Prerequisites:",
        "  cmake must be installed on your system (tested with ver 2.8.12.1)",
        "  the WT source must be in ./target/ios/wt/wt-3.3.4"
        ]);
end

operations = []

i = 0
while i < ARGV.length
    arg = ARGV[i]
    i += 1
    if (arg == '-help') then
        show_usage()
        exit(1)
    end
    if (arg == '-noredirect') then
        $REDIRECT = false
        next
    end
    operations.push(arg)
end

if operations.length == 0 then
    show_usage()
    exit 1
end


builder = Builder.new
preflight = PreFlight.new
preflight.check_directories

operations.each do |operation|
    if contains(operation, 'configtest') or contains(operation, 'all')
        preflight.check_xcode("6.4")
        preflight.check_sdk("iphoneos8.4")
    end

    operation = operation.downcase

    builder.build_all if operation == 'all'
    builder.clean_all if operation == 'clean'
    builder.rebuild_boost if operation == 'boost'
    builder.clean_wt if operation == 'wt'
    builder.rebuild_wt if operation == 'wt'
    #    builder.copy_wt_files if operation == 'wt'
    #    builder.copy_wt_files if operation == 'wtcopy'
    builder.clean_interspec if operation == 'interspec'
    builder.build_interspec if operation == 'interspec'
    builder.create_xcode_directories if operation == 'xcodedirs'
end
