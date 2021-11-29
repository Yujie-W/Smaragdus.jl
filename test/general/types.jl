@testset "General --- type constructors" begin
    for FT in [Float32, Float64]
        # create SoilAir struct
        sa = SoilAir{FT}(collect(FT,-2:0.2:0), collect(FT,0:0.2:4));
        @test true;
    end
end
