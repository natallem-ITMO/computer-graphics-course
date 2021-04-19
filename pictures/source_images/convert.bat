echo off
set arg1=%1
set arg2=%2
set arg3=%3
magick convert %arg1%.%arg2% "B:\Projects\GitProjects\Graphics\pictures\input_pictures\%arg1%.%arg3%"
echo on