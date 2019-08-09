using NonlinearEigenproblems;
using Test

@testset "dtn_dimer" begin
    # This requires that the FEM-files are downloaded
    nep=nep_gallery("dtn_dimer","/home/jarl/archive/dtn_umea_collab_nosync_julia/dimer_TM_ref2_p10_g0_a1",40);


    n=size(nep,1);
    # Reference eigvals computed with MATLAB code
    λv=[4.370036082655539 - 1.526561686254116im
        2.165520379441414 - 0.537312290190438im
        3.958068685793646 - 0.528959556453958im
        4.397326877370154 - 0.671081974763288im
        1.098616611079115 - 1.005745095689971im]

    λ=λv[end];


    # Solve it:
    (λs,vs)=quasinewton(nep,λ=1.099+1.006im,v=ones(n),logger=displaylevel,tol=1e-5)
    (λss,vss)=quasinewton(nep,λ=λs,v=vs,logger=displaylevel)

    @test minimum(abs.(λss.-λv))<1e-10 # gives 1.4512316607747323e-12

end
