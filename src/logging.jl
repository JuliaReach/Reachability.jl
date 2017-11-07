import Base: info, warn, toc
import Memento: debug

export info, warn, debug, tocc, configure_logger

global G_LOGGER = nothing

const DEFAULT_LOG_LEVEL = "warn"

@inline function info(msg::String)
    Memento.info(G_LOGGER, msg)
end

@inline function warn(msg::String)
    Memento.warn(G_LOGGER, msg)
end

@inline function debug(msg::String)
    Memento.debug(G_LOGGER, msg)
end

@inline function tocc(func::Function=info)
    func("elapsed time: $(toq()) seconds")
end

function configure_logger(level::Union{String, Int, Void}=DEFAULT_LOG_LEVEL)
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
    return Memento.config(level_string; fmt="[{level}] {msg}")
end
