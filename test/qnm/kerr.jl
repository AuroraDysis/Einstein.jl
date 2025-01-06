@testitem "qnm_kerrchebnep" begin
    using GRSuite, Test

    a = 0.7
    s = 2
    l = 2
    m = 2
    n = 80
    ω = 0.532600243551018 - 0.08079287315500905im
    tol = 1e-12

    cache = qnm_kerr_chebnep_cache(Float64, a, s, m, n)
    δ = qnm_kerr_chebnep_step!(cache, ω, l)

    @test abs(δ) < tol
end

@testitem "qnm_kerr" begin
    for TR in [Float64, BigFloat]
        a = parse(TR, "0.7")
        s = 2
        l = 2
        m = 2
        n = 0
        ω_re = parse(
            TR,
            "0.5326002435510180172081539794333463691995588941440439727628150162718570815930126401361418412332984881640541691941385982263104777733734855994767002752016226141820059433197419654377634139360200659490595368156399402599480889331192121863934312384706603135802905665249771764508297527804314693809975204195436081342474460458243417402921820003003496127253443889928798711683581822506847471481697238407948180554",
        )
        ω_im = parse(
            TR,
            "-0.08079287315500766015888115743900604319267896786014540887670013112649052854133103823795873692396073311299271045640874574464049500229173565861003095831673606400431868873221296509453440459055802446845133365313994374953348723521731852638597327067090339342934114665321022545244251053933463266555851966283212917828336870308396316090880733935787815211441638976043728595827022167896909165405672270564221308464",
        )
        ω = ω_re + ω_im * 1im
        A_re = parse(
            TR,
            "-1.0968330181260638646138712824950021525474649293724341226983465304343882430034856404501199519985417476267452199006946882255621463106887406919663249530149412577038302364739112497269205450888600677654417485951426219620991417879814180448745617489537225771920909555519463575794333725993529131342273595849310636681857148895076447889072622660185988265686593762188957862651797349582433687853401936694368834397900468180844141939507026",
        )
        A_im = parse(
            TR,
            "0.18315825271747034394648312674902077839702438802934906156891059819213745012812016631885215847554742867115659589113159212316248086425302135432085755384361760069916616285701635400467880482516485463646453921134473840339034357437627872414030755530030916707848877949515409238351034872189456032412626109394579740107729060210654091702942338338955851934908031025555585224521892558344818865147030116930534874228299066119736732483983",
        )
        A = A_re + A_im * 1im
        C = [
            0.9974661860984922 + 0.0im,
            -0.07022978582035719 + 0.010660390441034171im,
            0.003739828906170969 + -0.001151576549206862im,
            -0.00015181121425804852 + 7.292590142440934e-05im,
            4.944215061024621e-06 + -3.360651603310845e-06im,
            -1.3166989201504968e-07 + 1.2181467892545735e-07im,
            2.914709438277216e-09 + -3.644525568057317e-09im,
            -5.3573339010985764e-11 + 9.2550793070661e-11im,
            8.018496860593015e-13 + -2.038388358866514e-12im,
            -9.030840491596909e-15 + 3.950667194618473e-14im,
            5.097536947576771e-17 + -6.818514653125067e-16im,
            7.999510397309874e-19 + 1.0567762252271612e-17im,
            -3.400854559163453e-20 + -1.4806805794773245e-19im,
            7.4367915608113155e-22 + 1.8836431789635946e-21im,
            -1.26728822621701e-23 + -2.1813384779695933e-23im,
            1.846603234484284e-25 + 2.29936832495328e-25im,
            -2.389778913779243e-27 + -2.1989586577769056e-27im,
            2.801902062759508e-29 + 1.8901843833280047e-29im,
            -3.0124822517651177e-31 + -1.4296878207850558e-31im,
        ]
        l_max = 20
        ω_pert = ω + rand(Complex{TR}) / 10
        params = QNMKerrParams{TR}(; a=a, s=s, l=l, m=m, n=n, ω_guess=ω_pert, l_max=l_max)

        ωsol = qnm_kerr(params)

        tol = typetol(TR)
        min_tol = 10^(-20)
        tol = max(100 * tol, min_tol)
        @test abs(ωsol - ω) < tol
    end
end