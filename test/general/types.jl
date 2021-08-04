@testset "General --- type constructors" begin
    for FT in [Float32, Float64]
        _lbio = Emerald.LeafBio{FT}();
        _leaf = Leaf{FT}();
        _wls  = WaveLengthSet{FT}();
        @test true;
    end
end
