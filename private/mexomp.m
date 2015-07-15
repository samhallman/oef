function mexomp( cppFile )
% mexomp( cppFile )
% Compile mex file with OpenMP opts set

cmd = ['mex %s ''-DUSEOMP'' CXXFLAGS="\\$CXXFLAGS ' ...
       '-fopenmp" LDFLAGS="\\$LDFLAGS -fopenmp"'];
eval( sprintf(cmd,cppFile) );
