# ADAPTED BY ROMIN TOMASETTI - romin.tomasetti at gmail.com- 5th APRIL 2018
# From https://github.com/maidsafe-archive/MaidSafe/blob/master/cmake_modules/add_boost.cmake

set(BoostVersion 1.66.0)
if(MSVC)
	set(BoostZipExtension ".zip")
	set(BoostSHA256       e1c55ebb00886c1a96528e4024be98a38b815115f62ecfe878fcf587ba715aad)
	set(BoostDownloadLink "https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.zip")
else()
	set(BoostZipExtension ".tar.bz2")
	set(BoostSHA256       5721818253e6a0989583192f96782c4a98eb6204965316df9f5ad75819225ca9)
	set(BoostDownloadLink "https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2")
endif()

set(BoostCacheDir "${PROJECT_SOURCE_DIR}")

################################################################################
# Create build folder name :
set(BoostFolderName "${BoostVersion}")
string(REPLACE "." "_" BoostFolderName ${BoostFolderName})
set(BoostFolderName boost_${BoostFolderName})
message(STATUS "BoostFolderName set to ${BoostFolderName}")
################################################################################
# Set up the full path to the source directory
set(BoostSourceDir "${BoostFolderName}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}")
if(HAVE_LIBC++)
  set(BoostSourceDir "${BoostSourceDir}_LibCXX")
endif()
if(HAVE_LIBC++ABI)
  set(BoostSourceDir "${BoostSourceDir}_LibCXXABI")
endif()
if(CMAKE_CL_64)
  set(BoostSourceDir "${BoostSourceDir}_Win64")
endif()
if(${ANDROID_BUILD})
  set(BoostSourceDir "${BoostFolderName}_Android_v${AndroidApiLevel}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}")
endif()
string(REPLACE "." "_" BoostSourceDir ${BoostSourceDir})
set(BoostSourceDir "${BoostCacheDir}/${BoostSourceDir}")
message(STATUS "BoostSourceDir set to ${BoostSourceDir}")
################################################################################
# Check the full path to the source directory is not too long for Windows.  File paths must be less
# than MAX_PATH which is 260.  The current longest relative path Boost tries to create is:
# Build\boost\bin.v2\libs\program_options\build\fd41f4c7d882e24faa6837508d6e5384\libboost_program_options-vc120-mt-gd-1_55.lib.rsp
# which along with a leading separator is 129 chars in length.  This gives a maximum path available
# for 'BoostSourceDir' as 130 chars.
if(WIN32)
  get_filename_component(BoostSourceDirName "${BoostSourceDir}" NAME)
  string(LENGTH "/${BoostSourceDirName}" BoostSourceDirNameLengthWithSeparator)
  math(EXPR AvailableLength 130-${BoostSourceDirNameLengthWithSeparator})
  string(LENGTH "${BoostSourceDir}" BoostSourceDirLength)
  if(BoostSourceDirLength GREATER "130")
    set(Msg "\n\nThe path to boost's source is too long to handle all the files which will ")
    set(Msg "${Msg}be created when boost is built.  To avoid this, set the CMake variable ")
    set(Msg "${Msg}USE_BOOST_CACHE to ON and set the variable BOOST_CACHE_DIR to a path ")
    set(Msg "${Msg}which is at most ${AvailableLength} characters long.  For example:\n")
    set(Msg "${Msg}  mkdir C:\\maidsafe_boost\n")
    set(Msg "${Msg}  cmake . -DUSE_BOOST_CACHE=ON -DBOOST_CACHE_DIR=C:\\maidsafe_boost\n\n")
    message(FATAL_ERROR "${Msg}")
  endif()
endif()
################################################################################
# Download boost if required
set(ZipFilePath "${BoostCacheDir}/build/${BoostFolderName}${BoostZipExtension}")
if(NOT EXISTS ${ZipFilePath})
  message(STATUS "Downloading boost ${BoostVersion} to ${ZipFilePath}")
endif()
file(DOWNLOAD ${BoostDownloadLink}
     ${ZipFilePath}
     STATUS Status
     SHOW_PROGRESS
     EXPECTED_HASH SHA256=${BoostSHA256}
     )
################################################################################
# Extract boost if required
string(FIND "${Status}" "returning early" Found)
if(Found LESS "0" OR NOT IS_DIRECTORY "${BoostSourceDir}")
  set(BoostExtractFolder "${BoostCacheDir}/boost_unzip")
  file(REMOVE_RECURSE ${BoostExtractFolder})
  file(MAKE_DIRECTORY ${BoostExtractFolder})
  file(COPY ${ZipFilePath} DESTINATION ${BoostExtractFolder})
  message(STATUS "Extracting boost ${BoostVersion} to ${BoostExtractFolder}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${BoostFolderName}.tar.bz2
                  WORKING_DIRECTORY ${BoostExtractFolder}
                  RESULT_VARIABLE Result
                  )
  if(NOT Result EQUAL "0")
    message(FATAL_ERROR "Failed extracting boost ${BoostVersion} to ${BoostExtractFolder}")
  endif()
  file(REMOVE ${BoostExtractFolder}/${BoostFolderName}.tar.bz2)

  # Get the path to the extracted folder
  file(GLOB ExtractedDir "${BoostExtractFolder}/*")
  list(LENGTH ExtractedDir n)
  if(NOT n EQUAL "1" OR NOT IS_DIRECTORY ${ExtractedDir})
    message(FATAL_ERROR "Failed extracting boost ${BoostVersion} to ${BoostExtractFolder}")
  endif()
  file(RENAME ${ExtractedDir} ${BoostSourceDir})
  file(REMOVE_RECURSE ${BoostExtractFolder})
endif()
################################################################################
include(ProcessorCount)
ProcessorCount(N)
message("Your computer has ${N} cores.")
if(NOT IS_DIRECTORY "${BoostSourceDir}/bin.v2")
set(BOOTSTRP_COMMAND1 "--prefix=${BoostSourceDir}")
set(BOOTSTRP_COMMAND2 "--with-libraries=filesystem,regex,system")
execute_process(COMMAND ${BoostSourceDir}/bootstrap.sh ${BOOTSTRP_COMMAND1}
		WORKING_DIRECTORY ${BoostSourceDir})
execute_process(COMMAND ${BoostSourceDir}/bootstrap.sh ${BOOTSTRP_COMMAND2}
		WORKING_DIRECTORY ${BoostSourceDir})
execute_process(COMMAND ${BoostSourceDir}/b2 -j ${N}
		WORKING_DIRECTORY ${BoostSourceDir})
else()
	message(STATUS "Boost already exists and is already compiled.")
endif()


