import Base: info, warn, toc
import Memento: debug

export info, warn, debug, tocc,
       configure_logger, add_file_logger

global LOGGER = Memento.getlogger(current_module())

const DEFAULT_LOG_LEVEL = "warn"

"""
    info(msg)

Prints a message on info level.

### Input

- `msg` - message string
"""
@inline function info(msg::String)
    Memento.info(LOGGER, msg)
end

"""
    warn(msg)

Prints a message on warn level.

### Input

- `msg` - message string
"""
@inline function warn(msg::String)
    Memento.warn(LOGGER, msg)
end

"""
    debug(msg)

Prints a message on debug level.

### Input

- `msg` - message string
"""
@inline function debug(msg::String)
    Memento.debug(LOGGER, msg)
end

"""
    tocc(func)

Prints the output of `toc()` using the given log function (default: `info`).

### Input

- `func` - (optional, default: `info`) log function
"""
@inline function tocc(func::Function=info)
    func("elapsed time: $(toq()) seconds")
end

"""
    configure_logger(level)

Configures the global log level. If no log level is passed, we use the log level
that is defined by the constant `DEFAULT_LOG_LEVEL`.

### Input

- `level` - (optional) the log level; can be either an integer between 0 and 2
            or a string that is supported by the
            [Memento.jl](https://invenia.github.io/Memento.jl/latest/man/intro.html#Logging-levels-1)
            package.
"""
function configure_logger(level::Union{String, Int, Nothing}=DEFAULT_LOG_LEVEL)
    if level isa String
        level_string = level
    elseif level isa Int
        if level == 0
            level_string = "warn"
        elseif level == 1
            level_string = "info"
        elseif level == 2
            level_string = "debug"
        else
            error("Illegal verbosity input $level.")
        end
    else
        error("Illegal verbosity input $level.")
    end
    return Memento.config!(level_string; fmt="[{level}] {msg}")
end

"""
    add_file_logger(filename)

Sets up an additional logger to a file.

### Input

- `filename` -- (optional, default: `tempname()`) the log file name
"""
function add_file_logger(filename::String=tempname())::Nothing
    formatter = Memento.DefaultFormatter("[{date}|{level}] {msg}")
    handler = Memento.DefaultHandler(filename, formatter)
    push!(LOGGER, handler)
    return nothing
end
