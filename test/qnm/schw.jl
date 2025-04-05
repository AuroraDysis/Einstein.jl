@testitem "qnm_schw_expansion_dolan_ottewill" begin
    @testset "Gravitational Modes" begin
        for s in [-2, 2]
            @testset "s = $s" begin
                let l = 2, n = 0
                    ω = 0.3736716844180419 - 0.0889623156889357im
                    ϖ = qnm_schw_expansion_dolan_ottewill(Float64, s, l, n)
                    @test ϖ isa ComplexF64
                    @test abs(ϖ - ω) / abs(ω) < 0.01
                end

                let l = 2, n = 1
                    ω = 0.3467109968791634 - 0.27391487529123487im
                    ϖ = qnm_schw_expansion_dolan_ottewill(Float64, s, l, n)
                    @test ϖ isa ComplexF64
                    @test abs(ϖ - ω) / abs(ω) < 0.02
                end
            end
        end
    end
end
