import Base: info, toc

export info, toc

global G_LOG = false

function info(msg::String)
    if G_LOG
        Base.println(msg)
    end
end

function toc()
    if G_LOG
        Base.toc()
    end
end