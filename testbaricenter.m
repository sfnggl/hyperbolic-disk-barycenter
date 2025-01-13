clear all;

addpath(genpath('./src'))

opts1 = struct("search_type", "sd", "field", "real");
opts2 = struct("search_type", "bb", "field", "real");

%switch(input("provide n of testcase\n"))
switch(1)
  case 1
    poincarediskbaricenter(opts1)
  case 2
    poincarediskbaricenter(opts2)
  otherwise
    return
end
