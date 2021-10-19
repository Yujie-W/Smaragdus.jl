@testset "General --- type constructors" begin
    for FT in [Float32, Float64]
        # wavelength dataset
        _wls = WaveLengthSet{FT}();
        _wls = WaveLengthSet{FT}(collect(FT,400:5:2500));
        _wls = WaveLengthSet{FT}(collect(FT,400:5:2500); opti=Emerald.OPTI_2017);
        @test true;

        # create leaf biophysics
        _lbio = Emerald.LeafBiophysics{FT}();
        _lbio = Emerald.LeafBiophysics{FT}(collect(FT,400:5:2400.1));
        _lbio = Emerald.LeafBiophysics{FT}(collect(FT,400:5:2400.1); opti=Emerald.OPTI_2021);
        _lbio = Emerald.LeafBiophysics{FT}(WaveLengthSet{FT}());
        _lbio = Emerald.LeafBiophysics{FT}(WaveLengthSet{FT}(); opti=Emerald.OPTI_2021);
        @test true;

        # create leaf biophysics : warning expected here
        @info "A warning is expected here...";
        _lbio = Emerald.LeafBiophysics{FT}(collect(FT,400:5:2410.1));
        @test true;

        # create leaf
        _leaf = Leaf{FT}();
        _leaf = Leaf{FT}(collect(FT,400:5:2400.1));
        _leaf = Leaf{FT}(collect(FT,400:5:2400.1); opti=Emerald.OPTI_2021);
        _leaf = Leaf{FT}(WaveLengthSet{FT}());
        _leaf = Leaf{FT}(WaveLengthSet{FT}(); opti=Emerald.OPTI_2021);
        @test true;

        # create van Genuchten soil
        _vg = VanGenuchten{FT}("Loam");
        _vg = VanGenuchten{FT}("Silt");
        _vg = VanGenuchten{FT}("Test", FT(100), FT(2), FT(0.5), FT(0.1));
        @test true;

        # create van Genuchten soil : warning expected here
        @info "A warning is expected here...";
        _vg = VanGenuchten{FT}("Unsupported_Soil");
        @test true;

        # create Brooks Corey soil
        bc = BrooksCorey{FT}("Test", FT(5), FT(2), FT(0.5), FT(0.1));
        bc = BrooksCorey{FT}(VanGenuchten{FT}("Loam"));
        @test true;

        # create SoilAir struct
        sa = SoilAir{FT}(collect(FT,-2:0.2:0), collect(FT,0:0.2:4));
        @test true;
    end
end
