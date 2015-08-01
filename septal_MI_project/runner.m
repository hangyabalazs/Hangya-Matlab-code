fix(clock)
try
    b_entryrun4
    fix(clock)
catch
    fix(clock)
    rethrow(lasterror)
end