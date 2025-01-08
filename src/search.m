function [xs, ds, steps] = search(f,g,x0,opts)

  %% define possible search types
  possible_search_types = [{"sd" "Sd" "SD"};{"bb" "Bb" "BB"}];

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
  switch (opts.search_type)
    case possible_search_types(1,:)
      search_type = "sd";
    case possible_search_types(2,:)
      search_type = "bb";
    otherwise
      disp ("possible values are: ")
      disp (possible_search_types)
      error ("invalid search type")
  end

  %% begin search and output results
  tic;
  switch (search_type)
    case {"sd"}
      [xs, ds, steps] = sd(f,g,x0,max_iter,tol)
    case {"bb"}
      [xs, ds, steps] = bb(f,g,x0,max_iter,tol)
    otherwise
      error ("invalid search type")
  end
  toc
end

