function [xs, ds, steps] = search(f,g,x0,search_type)

  %% check standard arguments
  if ~exist("x0")
    error("no starting point provided.")
  elseif ~exist("f")
    error("no function to evaluate provided.")
  elseif ~exist("g")
    error("no gradient to evaluate provided.")
  end

  %% begin search and output results
  switch (search_type)
    case {"sd"}
      tic;
      [xs, ds, steps] = sd(f,g,x0);
      toc
    case {"bb"}
      tic;
      [xs, ds, steps] = bb(f,g,x0);
      toc
    otherwise
      error("Invalid search type.\n\
Possible search types are: sd, bb\n\n")
      exit()
  end
end

