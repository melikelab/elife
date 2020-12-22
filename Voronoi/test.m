 % load the C++ library needed to read .nd2 quickly
            dir = 'W:\AFIB_SOFTWARE\Matlab_lib\io\nd2lib';
            libname = fullfile(dir, dllfile);
            hfile = fullfile(dir, 'nd2ReadSDK.h');
            loadlibrary(libname, hfile);