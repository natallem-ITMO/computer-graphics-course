echo off
set arg1=%1
set arg2=%2
magick convert %arg1%.%arg2% %arg1%_.pgm
echo on