# disc cross section profiles from Wiktoria's work
# add corona location
# fix the aspect ratio to show that the disc is actually quite thin

using Gradus, Plots, Printf
using LaTeXStrings
r_start=-50
r_end=50
n=r_end-r_start
r=LinRange(r_start,r_end,n)

function Gradus.cross_section(d, ρ)
    y=Gradus.isco(m)
    if ρ < y
        return -one(typeof(ρ))
    end
    H = (3 / 2) * inv(η) * (d.Ṁ / d.Ṁedd) * (1 - sqrt(y / ρ))
end

function all_height(disk)
    height=zeros(r_end)
    height_b=zeros(r_end)
    for i in 1:r_end
        height[i]=cross_section(disk,i)
        height_b[i]=0-cross_section(disk,i)
    end
    a=reverse(height)
    b=reverse(height_b)
    top=vcat(a,height)
    bottom=vcat(b,height_b) 
    return top, bottom
end

function disk_height(m)
    d1=ShakuraSunyaev(m,eddington_ratio=0.1)
    d2=ShakuraSunyaev(m,eddington_ratio=0.2)
    d3=ShakuraSunyaev(m,eddington_ratio=0.3)

    d1_t,d1_b=all_height(d1)
    d2_t,d2_b=all_height(d2)
    d3_t,d3_b=all_height(d3)
    return d1_t,d1_b,d2_t,d2_b,d3_t,d3_b
end

function plot_data(dx_at,dx_ab,dx_bt,dx_bb,dx_ct,dx_cb)
    p=plot(r,[dx_ct],fillrange=dx_cb,label=L"\dot{M} / \dot{M}_{\textrm{Edd}} = 0.3", color = :darkorchid1, fillcolor = :darkorchid1, xlabel=L"\textrm{Radius~} (GM/c^2)",ylabel=L"\textrm{Height~} (GM/c^2)", aspect_ratio=1, xrange=(-25,25), yrange=(-10,10), size=(600,300))
    plot!(r,dx_bt,fillrange=dx_bb,label=L"\dot{M} / \dot{M}_\textrm{Edd} = 0.2", color = :hotpink2, fillcolor = :hotpink2)
    plot!(r,dx_at,fillrange=dx_ab, label=L"\dot{M} / \dot{M}_\textrm{Edd} = 0.1", color = :tan1, fillcolor = :tan1)
    x=zeros(100)
    plot!(r,x,label=L"\dot{M} / \dot{M}_\textrm{Edd} = 0.0",color=:white)
    return p
end  
#a=0.0 α13=0.35
# m1_a=JohannsenMetric(M=1.0,a=0.0,α13=0.0,ϵ3=0.0)
m1_a=KerrMetric(1.0,0.0)
d1_at,d1_ab,d1_bt,d1_bb,d1_ct,d1_cb=disk_height(m1_a)
plot_data(d1_at,d1_ab,d1_bt,d1_bb,d1_ct,d1_cb)
plot_horizon!(m1_a, fill=true, color = :black, label = L"\textrm{BH}")
plot!([0], [8], marker = :star, markersize = 12, color = :red, linecolor = :white, label = L"\textrm{Corona}")
savefig("presentation/disc_profile_a_0.svg")
