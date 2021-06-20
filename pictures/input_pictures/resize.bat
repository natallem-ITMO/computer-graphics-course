echo off
set name=%1
set w=%2
set h=%3
magick convert  %name% -gamma 0.454545 -filter Box -resize %w%x%h%! -gamma 2.2 "B:\Projects\GitProjects\Graphics\pictures\output_pictures\magick\%name%"
echo on
