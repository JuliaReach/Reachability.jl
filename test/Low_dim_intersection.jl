using LazySets, Optim
import LazySets.Approximations:decompose, project, overapproximate

n0 = 16
z = zeros(n0)

X0 = Hyperrectangle(low=[0.;0.],
                    high=[1.; 1.]);

proj_X0 = decompose(X0, blocks=[1, dim(X0)-1])
G = HPolyhedron([HalfSpace([1.0; 0.], 0.5)])
proj_G = project(G, [1], LinearMap)

proj_inter = proj_G ∩ proj_X0.array[1]
proj_X0.array[1] = overapproximate(proj_inter, Hyperrectangle)

plot(project(proj_X0,[1,2], Hyperrectangle))
