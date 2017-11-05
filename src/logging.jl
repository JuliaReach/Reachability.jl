import Base: info, toc

export info, tocc

global G_LOG = false

function info(msg::String)
    if G_LOG
        println(msg)
    end
end

function tocc()
    if G_LOG
        toc()
    end
end