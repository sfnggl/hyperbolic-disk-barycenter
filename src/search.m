function [xs, ds, steps] = search(f,g,x0,opts)

  %% check standard arguments
  if ~exist("x0")
    error("no starting point provided.")
  elseif ~exist("f")
    error("no function to evaluate provided.")
  elseif ~exist("g")
    error("no gradient to evaluate provided.")
  elseif ~exist("opts")
    search_type = "sd";
    % initialize maximum iteration and tolerance if not given
  end

  %% check for search type specified by user
  search_type = opts.search_type;

  %% begin search and output results
  switch (search_type)
    case {"sd"}
      tic;
      [xs, ds, steps] = sd(f,g,x0)
      toc
    case {"bb"}
      tic;
      [xs, ds, steps] = bb(f,g,x0)
      toc
    otherwise
      disp(...
	    "Invalid search type.
	    Possible search types are: sd, bb")
      exit()
  end
end

