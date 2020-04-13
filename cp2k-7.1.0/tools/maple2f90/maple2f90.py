#!/bin/env python
import re,sys
"""Script to help the conversion of maple generated fortran code to f90"""

renamedVar=["cg", "cg0", "cg1", "cg10", "cg11", "cg12", "cg13", "cg14", "cg15", "cg16", "cg17", "cg18", "cg19", "cg2", "cg20", "cg21", "cg22", "cg23", "cg24", "cg25", "cg26", "cg27", "cg28", "cg29", "cg3", "cg30", "cg31", "cg32", "cg33", "cg34", "cg35", "cg36", "cg37", "cg38", "cg39", "cg4", "cg40", "cg41", "cg42", "cg43", "cg44", "cg45", "cg46", "cg47", "cg48", "cg49", "cg5", "cg50", "cg51", "cg52", "cg53", "cg54", "cg55", "cg56", "cg57", "cg58", "cg59", "cg6", "cg60", "cg61", "cg62", "cg63", "cg64", "cg65", "cg66", "cg67", "cg68", "cg69", "cg7", "cg70", "cg71", "cg72", "cg73", "cg74", "cg75", "cg76", "cg77", "cg78", "cg79", "cg8", "cg80", "cg81", "cg82", "cg83", "cg84", "cg85", "cg86", "cg87", "cg88", "cg89", "cg9", "cg90", "cg91", "cg92", "cg93", "cg94"]

origNames=["norm_drho", "norm_drhoa", "norm_drhob", "chirhoa", "epsilon_c_unifrhoarhob", "e_c_u_01rhob", "epsilon_c_unif1rhob", "alpha_c1rhoa", "k_frhoarhob", "chirhobrhob", "phi1rhob", "phirhoa", "kf_b", "kf_brhob", "Fx_b", "e_c_u_01rhoa", "frhoarhob", "chirhoarhob", "chirhob", "Fx_a", "ex_unif_a1rhoa", "e_c_u_1rhoa", "epsilon_cGGArhob", "epsilon_c_unifrhob", "k_frhoa", "phirhoarhoa", "trhobrhob", "ex_unif_b", "Arhoarhob", "k_s1rhoa", "gamma_var", "k_frhob", "epsilon_cGGA", "phirhob", "Fx_bnorm_drhob", "s_anorm_drhoa", "epsilon_cGGArhoa", "Arhobrhob", "Fx_brhob", "k_s1rhob", "s_arhoa", "rsrhoarhoa", "epsilon_c_unif", "kf_arhoarhoa", "Hnorm_drho", "trhoarhoa", "Fx_anorm_drhoa", "alpha_crhob", "e_c_u_0rhob", "s_a1rhoa", "k_s", "kf_a", "k_srhoa", "rsrhoarhob", "rsrhobrhob", "s_a", "epsilon_c_unifrhoa", "e_c_u_0rhoa", "phi1rhoa", "trhoarhob", "alpha_c", "frhoarhoa", "trhobnorm_drho", "Fx_a1rhoa", "e_c_u_0rhobrhob", "Arhoarhoa", "phirhoarhob", "ex_unif_brhob", "kf_arhoa", "frhobrhob", "phirhobrhob", "Fx_arhoa", "e_c_u_0rhoarhob", "s_bnorm_drhob", "epsilon_c_unifrhoarhoa", "epsilon_c_unif1rhoa", "e_c_u_0rhoarhoa", "ex_unif_arhoa", "epsilon_c_unifrhobrhob", "s_brhob", "k_srhob", "ex_unif_a", "chirhoarhoa", "e_c_u_1rhob", "s_b", "alpha_c1rhob", "tnorm_drho", "k_frhoarhoa", "s_b1rhob", "kf_brhobrhob", "e_c_u_0", "ex_unif_b1rhob", "trhoanorm_drho", "Fx_b1rhob", "alpha_crhoa", "r_eqs_lsd4"]

code="""
      doubleprecision function cg94 (rhoa, rhob, cg, cg0, cg1)

        my_rhoa = rhoa
        my_rhob = rhob
        my_norm_drho = cg
        my_norm_drhoa = cg0
        my_norm_drhob = cg1
        rho = rhoa + rhob
        t1 = rhoa - rhob
        t2 = 0.1D1 / rho
        chi = t1 * t2
        t3 = 3 ** (0.1D1 / 0.3D1)
        t4 = 4 ** (0.1D1 / 0.3D1)
        t5 = t4 ** 2
        t6 = t3 * t5
        t7 = 0.1D1 / 0.3141592654D1
        t8 = t7 * t2
        t9 = t8 ** (0.1D1 / 0.3D1)
        rs = dble(t6) * t9 / 0.4D1
        t12 = 0.1D1 + 0.21370D0 * rs
        t13 = sqrt(rs)
        t16 = t13 * rs
        t18 = rs ** 0.20D1
        t20 = 0.75957D1 * t13 + 0.35876D1 * rs + 0.16382D1 * t16 + 0.492
     #94D0 * t18
        t23 = 0.1D1 + 0.1608182432D2 / t20
        t24 = log(t23)
        cg9 = -0.62182D-1 * t12 * t24
        t28 = 0.1D1 + 0.20548D0 * rs
        t33 = 0.141189D2 * t13 + 0.61977D1 * rs + 0.33662D1 * t16 + 0.62
     #517D0 * t18
        t36 = 0.1D1 + 0.3216468318D2 / t33
        t37 = log(t36)
        t41 = 0.1D1 + 0.11125D0 * rs
        t46 = 0.10357D2 * t13 + 0.36231D1 * rs + 0.88026D0 * t16 + 0.496
     #71D0 * t18
        t49 = 0.1D1 + 0.2960857464D1 / t46
        t50 = log(t49)
        cg62 = 0.33774D0 * t41 * t50
        t52 = 2 ** (0.1D1 / 0.3D1)
        t55 = 1 / (2 * t52 - 2)
        t56 = 0.1D1 + chi
        t57 = t56 ** (0.1D1 / 0.3D1)
        t58 = t57 * t56
        t59 = 0.1D1 - chi
        t60 = t59 ** (0.1D1 / 0.3D1)
        t61 = t60 * t59
        f = (t58 + t61 - 0.2D1) * dble(t55)
        t63 = cg62 * f
        t64 = 0.9D1 / 0.8D1 / dble(t55)
        t65 = chi ** 2
        t66 = t65 ** 2
        t68 = t64 * (0.1D1 - t66)
        t70 = -0.31090D-1 * t28 * t37 - cg9
        t71 = t70 * f
        cg46 = cg9 + t63 * t68 + t71 * t66
        t73 = log(0.2D1)
        t75 = 0.3141592654D1 ** 2
        t76 = 0.1D1 / t75
        cg35 = (0.1D1 - t73) * t76
        t77 = t57 ** 2
        t78 = t60 ** 2
        phi = t77 / 0.2D1 + t78 / 0.2D1
        t80 = t75 * rho
        t81 = t80 ** (0.1D1 / 0.3D1)
        t82 = dble(t3) * t81 * t7
        t83 = sqrt(t82)
        cg53 = 0.2D1 * t83
        t84 = 0.1D1 / phi
        t85 = cg * t84
        t86 = 0.1D1 / cg53
        t87 = t86 * t2
        t = t85 * t87 / 0.2D1
        t89 = 0.1D1 / cg35
        t90 = cg46 * t89
        t91 = phi ** 2
        t92 = t91 * phi
        t93 = 0.1D1 / t92
        t95 = exp(-t90 * t93)
        t96 = t95 - 0.1D1
        A = 0.66725D-1 * t89 / t96
        t99 = cg35 * t92
        t100 = t ** 2
        t101 = t89 * t100
        t102 = A * t100
        t103 = 0.1D1 + t102
        t104 = A ** 2
        t105 = t100 ** 2
        t107 = 0.1D1 + t102 + t104 * t105
        t108 = 0.1D1 / t107
        t109 = t103 * t108
        t112 = 0.1D1 + 0.66725D-1 * t101 * t109
        t113 = log(t112)
        cg37 = cg46 + t99 * t113
        mu = 0.2224166667D-1 * t75
        t114 = t3 * t52
        t115 = t75 * rhoa
        t116 = t115 ** (0.1D1 / 0.3D1)
        cg54 = dble(t114) * t116
        cg81 = -0.3D1 / 0.4D1 * t7 * cg54
        t119 = 0.1D1 / cg54
        t120 = cg0 * t119
        t121 = 0.1D1 / rhoa
        cg58 = t120 * t121 / 0.2D1
        t123 = cg58 ** 2
        t126 = 0.1D1 + 0.1243781095D1 * mu * t123
        cg25 = 0.1804D1 - 0.804D0 / t126
        t129 = rhoa * cg81
        t131 = t75 * rhob
        t132 = t131 ** (0.1D1 / 0.3D1)
        cg19 = dble(t114) * t132
        cg32 = -0.3D1 / 0.4D1 * t7 * cg19
        t135 = 0.1D1 / cg19
        t136 = cg1 * t135
        t137 = 0.1D1 / rhob
        cg84 = t136 * t137 / 0.2D1
        t139 = cg84 ** 2
        t142 = 0.1D1 + 0.1243781095D1 * mu * t139
        cg20 = 0.1804D1 - 0.804D0 / t142
        t145 = rhob * cg32
        exc = t129 * cg25 + t145 * cg20 + rho * cg37
        t149 = rho ** 2
        t150 = 0.1D1 / t149
        t151 = t1 * t150
        cg10 = t2 - t151
        t152 = t9 ** 2
        t154 = 0.1D1 / t152 * t7
        rsrhoa = -dble(t6) * t154 * t150 / 0.12D2
        t160 = t20 ** 2
        t161 = 0.1D1 / t160
        t162 = t12 * t161
        t163 = 0.1D1 / t13
        t164 = t163 * rsrhoa
        t167 = t13 * rsrhoa
        t169 = rs ** 0.10D1
        t170 = t169 * rsrhoa
        t172 = 0.3797850000D1 * t164 + 0.35876D1 * rsrhoa + 0.2457300000
     #D1 * t167 + 0.985880D0 * t170
        t173 = 0.1D1 / t23
        t174 = t172 * t173
        cg6 = -0.1328829340D-1 * rsrhoa * t24 + 0.9999999999D0 * t162 *
     #t174
        t179 = t33 ** 2
        t180 = 0.1D1 / t179
        t181 = t28 * t180
        t186 = 0.7059450000D1 * t164 + 0.61977D1 * rsrhoa + 0.5049300000
     #D1 * t167 + 0.1250340D1 * t170
        t187 = 0.1D1 / t36
        t188 = t186 * t187
        cg27 = -0.638837320D-2 * rsrhoa * t37 + 0.1000000000D1 * t181 *
     #t188
        t193 = t46 ** 2
        t194 = 0.1D1 / t193
        t195 = t41 * t194
        t200 = 0.5178500000D1 * t164 + 0.36231D1 * rsrhoa + 0.1320390000
     #D1 * t167 + 0.993420D0 * t170
        t201 = 0.1D1 / t49
        t202 = t200 * t201
        cg93 = 0.375735750D-1 * rsrhoa * t50 - 0.9999999999D0 * t195 * t
     #202
        frhoa = (0.4D1 / 0.3D1 * t57 * cg10 - 0.4D1 / 0.3D1 * t60 * cg10
     #) * dble(t55)
        t209 = cg93 * f
        t211 = cg62 * frhoa
        t213 = t65 * chi
        t214 = t64 * t213
        t215 = t214 * cg10
        t217 = 0.4D1 * t63 * t215
        t218 = cg27 - cg6
        t219 = t218 * f
        t221 = t70 * frhoa
        t223 = t213 * cg10
        t225 = 0.4D1 * t71 * t223
        cg59 = cg6 + t209 * t68 + t211 * t68 - t217 + t219 * t66 + t221
     #* t66 + t225
        t226 = 0.1D1 / t57
        t228 = 0.1D1 / t60
        cg18 = t226 * cg10 / 0.3D1 - t228 * cg10 / 0.3D1
        t231 = t81 ** 2
        cg3 = dble(t3) / t231 * t75 / 0.3D1
        t235 = 0.1D1 / t83
        cg55 = t235 * cg3 * t7
        t237 = 0.1D1 / t91
        t238 = cg * t237
        t241 = cg53 ** 2
        t242 = 0.1D1 / t241
        t243 = t242 * t2
        t246 = t86 * t150
        t247 = t85 * t246
        trhoa = -t238 * t87 * cg18 / 0.2D1 - t85 * t243 * cg55 / 0.2D1 -
     # t247 / 0.2D1
        t249 = t96 ** 2
        t251 = t89 / t249
        t252 = cg59 * t89
        t254 = t91 ** 2
        t255 = 0.1D1 / t254
        t256 = t255 * cg18
        t259 = -t252 * t93 + 0.3D1 * t90 * t256
        Arhoa = -0.66725D-1 * t251 * t259 * t95
        t263 = cg35 * t91
        t264 = t113 * cg18
        t267 = t89 * t
        t268 = t109 * trhoa
        t271 = Arhoa * t100
        t272 = A * t
        t274 = 0.2D1 * t272 * trhoa
        t275 = t271 + t274
        t276 = t275 * t108
        t279 = t107 ** 2
        t280 = 0.1D1 / t279
        t281 = t103 * t280
        t282 = A * t105
        t285 = t100 * t
        t286 = t104 * t285
        t289 = t271 + t274 + 0.2D1 * t282 * Arhoa + 0.4D1 * t286 * trhoa
        t293 = 0.133450D0 * t267 * t268 + 0.66725D-1 * t101 * t276 - 0.6
     #6725D-1 * t101 * t281 * t289
        t294 = 0.1D1 / t112
        t295 = t293 * t294
        cg40 = cg59 + 0.3D1 * t263 * t264 + t99 * t295
        t298 = t116 ** 2
        cg7 = dble(t114) / t298 * t75 / 0.3D1
        cg78 = -0.3D1 / 0.4D1 * t7 * cg7
        t304 = cg54 ** 2
        t305 = 0.1D1 / t304
        t306 = cg0 * t305
        t309 = rhoa ** 2
        t310 = 0.1D1 / t309
        cg44 = -t306 * t121 * cg7 / 0.2D1 - t120 * t310 / 0.2D1
        t313 = t126 ** 2
        t315 = 0.1D1 / t313 * mu
        cg72 = 0.2000000001D1 * t315 * cg58 * cg44
        t319 = rhoa * cg78
        exc_rhoa = cg81 * cg25 + t319 * cg25 + t129 * cg72 + cg37 + rho
     #* cg40
        cg24 = -t2 - t151
        rsrhob = rsrhoa
        t326 = t163 * rsrhob
        t329 = t13 * rsrhob
        t331 = t169 * rsrhob
        t333 = 0.3797850000D1 * t326 + 0.35876D1 * rsrhob + 0.2457300000
     #D1 * t329 + 0.985880D0 * t331
        t334 = t333 * t173
        cg51 = -0.1328829340D-1 * rsrhob * t24 + 0.9999999999D0 * t162 *
     # t334
        t343 = 0.7059450000D1 * t326 + 0.61977D1 * rsrhob + 0.5049300000
     #D1 * t329 + 0.1250340D1 * t331
        t344 = t343 * t187
        cg83 = -0.638837320D-2 * rsrhob * t37 + 0.1000000000D1 * t181 *
     #t344
        t353 = 0.5178500000D1 * t326 + 0.36231D1 * rsrhob + 0.1320390000
     #D1 * t329 + 0.993420D0 * t331
        t354 = t353 * t201
        cg50 = 0.375735750D-1 * rsrhob * t50 - 0.9999999999D0 * t195 * t
     #354
        frhob = (0.4D1 / 0.3D1 * t57 * cg24 - 0.4D1 / 0.3D1 * t60 * cg24
     #) * dble(t55)
        t361 = cg50 * f
        t363 = cg62 * frhob
        t365 = t214 * cg24
        t367 = 0.4D1 * t63 * t365
        t368 = cg83 - cg51
        t369 = t368 * f
        t371 = t70 * frhob
        t373 = t213 * cg24
        t375 = 0.4D1 * t71 * t373
        cg29 = cg51 + t361 * t68 + t363 * t68 - t367 + t369 * t66 + t371
     # * t66 + t375
        cg38 = t226 * cg24 / 0.3D1 - t228 * cg24 / 0.3D1
        cg36 = cg3
        cg80 = t235 * cg36 * t7
        trhob = -t238 * t87 * cg38 / 0.2D1 - t85 * t243 * cg80 / 0.2D1 -
     # t247 / 0.2D1
        t385 = cg29 * t89
        t387 = t255 * cg38
        t390 = -t385 * t93 + 0.3D1 * t90 * t387
        Arhob = -0.66725D-1 * t251 * t390 * t95
        t394 = t113 * cg38
        t397 = t109 * trhob
        t400 = Arhob * t100
        t402 = 0.2D1 * t272 * trhob
        t403 = t400 + t402
        t404 = t403 * t108
        t411 = t400 + t402 + 0.2D1 * t282 * Arhob + 0.4D1 * t286 * trhob
        t415 = 0.133450D0 * t267 * t397 + 0.66725D-1 * t101 * t404 - 0.6
     #6725D-1 * t101 * t281 * t411
        t416 = t415 * t294
        cg28 = cg29 + 0.3D1 * t263 * t394 + t99 * t416
        t419 = t132 ** 2
        cg2 = dble(t114) / t419 * t75 / 0.3D1
        cg69 = -0.3D1 / 0.4D1 * t7 * cg2
        t425 = cg19 ** 2
        t426 = 0.1D1 / t425
        t427 = cg1 * t426
        t430 = rhob ** 2
        t431 = 0.1D1 / t430
        cg8 = -t427 * t137 * cg2 / 0.2D1 - t136 * t431 / 0.2D1
        t434 = t142 ** 2
        t436 = 0.1D1 / t434 * mu
        cg42 = 0.2000000001D1 * t436 * cg84 * cg8
        t440 = rhob * cg69
        exc_rhob = cg32 * cg20 + t440 * cg20 + t145 * cg42 + cg37 + rho
     #* cg28
        t445 = t84 * t86
        cg86 = t445 * t2 / 0.2D1
        t450 = t89 * t285
        t451 = A * cg86
        t459 = 0.2D1 * t272 * cg86 + 0.4D1 * t286 * cg86
        t463 = 0.133450D0 * t267 * t109 * cg86 + 0.133450D0 * t450 * t45
     #1 * t108 - 0.66725D-1 * t101 * t281 * t459
        t464 = t463 * t294
        cg48 = t99 * t464
        exc_norm_drho = rho * cg48
        cg4 = t119 * t121 / 0.2D1
        cg5 = 0.2000000001D1 * t315 * cg58 * cg4
        exc_norm_drhoa = t129 * cg5
        cg74 = t135 * t137 / 0.2D1
        cg39 = 0.2000000001D1 * t436 * cg84 * cg74
        exc_norm_drhob = t145 * cg39
        t474 = 0.1D1 / t149 / rho
        t475 = t1 * t474
        cg82 = -0.2D1 * t150 + 0.2D1 * t475
        t480 = t149 ** 2
        cg45 = -dble(t6) / t152 / t8 * t76 / t480 / 0.18D2 + dble(t6) *
     #t154 * t474 / 0.6D1
        t490 = rsrhoa * t161
        t495 = t12 / t160 / t20
        t496 = t172 ** 2
        t500 = 0.1D1 / t16
        t501 = rsrhoa ** 2
        t502 = t500 * t501
        t504 = t163 * cg45
        t507 = t163 * t501
        t509 = t13 * cg45
        t512 = t169 * cg45
        t518 = t160 ** 2
        t520 = t12 / t518
        t521 = t23 ** 2
        t522 = 0.1D1 / t521
        cg77 = -0.1328829340D-1 * cg45 * t24 + 0.4274000000D0 * t490 * t
     #174 - 0.2000000000D1 * t495 * t496 * t173 + 0.9999999999D0 * t162
     #* (-0.1898925000D1 * t502 + 0.3797850000D1 * t504 + 0.35876D1 * cg
     #45 + 0.1228650000D1 * t507 + 0.2457300000D1 * t509 + 0.9858800D0 *
     # t501 + 0.985880D0 * t512) * t173 + 0.1608182432D2 * t520 * t496 *
     # t522
        cg21 = cg6
        t528 = rsrhoa * t180
        t533 = t28 / t179 / t33
        t534 = t186 ** 2
        t549 = t179 ** 2
        t551 = t28 / t549
        t552 = t36 ** 2
        t553 = 0.1D1 / t552
        t559 = rsrhoa * t194
        t564 = t41 / t193 / t46
        t565 = t200 ** 2
        t580 = t193 ** 2
        t582 = t41 / t580
        t583 = t49 ** 2
        t584 = 0.1D1 / t583
        cg14 = cg93
        t588 = 0.1D1 / t77
        t589 = cg10 ** 2
        t594 = 0.1D1 / t78
        cg63 = (0.4D1 / 0.9D1 * t588 * t589 + 0.4D1 / 0.3D1 * t57 * cg82
     # + 0.4D1 / 0.9D1 * t594 * t589 - 0.4D1 / 0.3D1 * t60 * cg82) * dbl
     #e(t55)
        f1rhoa = frhoa
        t612 = cg14 * f
        t615 = cg62 * f1rhoa
        t618 = t64 * t65
        t633 = cg27 - cg21
        t640 = t633 * f
        t643 = t70 * f1rhoa
        t652 = -0.4D1 * t63 * t214 * cg82 + (-0.638837320D-2 * cg45 * t3
     #7 + 0.4109600000D0 * t528 * t188 - 0.2000000000D1 * t533 * t534 *
     #t187 + 0.1000000000D1 * t181 * (-0.3529725000D1 * t502 + 0.7059450
     #000D1 * t504 + 0.61977D1 * cg45 + 0.2524650000D1 * t507 + 0.504930
     #0000D1 * t509 + 0.12503400D1 * t501 + 0.1250340D1 * t512) * t187 +
     # 0.3216468318D2 * t551 * t534 * t553 - cg77) * f * t66 + t218 * f1
     #rhoa * t66 + 0.4D1 * t219 * t223 + t633 * frhoa * t66 + t70 * cg63
     # * t66 + 0.4D1 * t221 * t223 + 0.4D1 * t640 * t223 + 0.4D1 * t643
     #* t223 + 0.12D2 * t71 * t65 * t589 + 0.4D1 * t71 * t213 * cg82
        cg75 = cg77 + (0.375735750D-1 * cg45 * t50 - 0.2225000000D0 * t5
     #59 * t202 + 0.2000000000D1 * t564 * t565 * t201 - 0.9999999999D0 *
     # t195 * (-0.2589250000D1 * t502 + 0.5178500000D1 * t504 + 0.36231D
     #1 * cg45 + 0.6601950000D0 * t507 + 0.1320390000D1 * t509 + 0.99342
     #00D0 * t501 + 0.993420D0 * t512) * t201 - 0.2960857464D1 * t582 *
     #t565 * t584) * f * t68 + cg93 * f1rhoa * t68 - 0.4D1 * t209 * t215
     # + cg14 * frhoa * t68 + cg62 * cg63 * t68 - 0.4D1 * t211 * t215 -
     #0.4D1 * t612 * t215 - 0.4D1 * t615 * t215 - 0.12D2 * t63 * t618 *
     #t589 + t652
        cg76 = cg21 + t612 * t68 + t615 * t68 - t217 + t640 * t66 + t643
     # * t66 + t225
        t657 = 0.1D1 / t58
        t662 = 0.1D1 / t61
        cg30 = -t657 * t589 / 0.9D1 + t226 * cg82 / 0.3D1 - t662 * t589
     #/ 0.9D1 - t228 * cg82 / 0.3D1
        cg60 = cg18
        t670 = t75 ** 2
        cg87 = -0.2D1 / 0.9D1 * dble(t3) / t231 / t80 * t670
        t674 = 0.1D1 / t83 / t82
        t675 = cg3 ** 2
        cg34 = cg55
        t682 = cg * t93 * t86
        t683 = t2 * cg18
        t686 = t238 * t242
        t692 = t238 * t246 * cg18 / 0.2D1
        t696 = t2 * cg55
        t702 = t85 / t241 / cg53
        t705 = t242 * t150
        t708 = t85 * t705 * cg55 / 0.2D1
        t719 = t85 * t86 * t474
        cg49 = t682 * t683 * cg60 + t686 * t683 * cg34 / 0.2D1 + t692 -
     #t238 * t87 * cg30 / 0.2D1 + t686 * t696 * cg60 / 0.2D1 + t702 * t6
     #96 * cg34 + t708 - t85 * t243 * (-t674 * t675 * t76 / 0.2D1 + t235
     # * cg87 * t7) / 0.2D1 + t238 * t246 * cg60 / 0.2D1 + t85 * t705 *
     #cg34 / 0.2D1 + t719
        t1rhoa = -t238 * t87 * cg60 / 0.2D1 - t85 * t243 * cg34 / 0.2D1
     #- t247 / 0.2D1
        t727 = t89 / t249 / t96
        t728 = t95 ** 2
        t729 = t259 * t728
        t730 = cg76 * t89
        t732 = t255 * cg60
        t735 = -t730 * t93 + 0.3D1 * t90 * t732
        t746 = 0.1D1 / t254 / phi
        t747 = t746 * cg18
        cg67 = 0.133450D0 * t727 * t729 * t735 - 0.66725D-1 * t251 * (-c
     #g75 * t89 * t93 + 0.3D1 * t252 * t732 + 0.3D1 * t730 * t256 - 0.12
     #D2 * t90 * t747 * cg60 + 0.3D1 * t90 * t255 * cg30) * t95 - 0.6672
     #5D-1 * t251 * t259 * t735 * t95
        A1rhoa = -0.66725D-1 * t251 * t735 * t95
        t765 = cg35 * phi
        t772 = A1rhoa * t100
        t774 = 0.2D1 * t272 * t1rhoa
        t775 = t772 + t774
        t776 = t775 * t108
        t783 = t772 + t774 + 0.2D1 * t282 * A1rhoa + 0.4D1 * t286 * t1rh
     #oa
        t787 = 0.133450D0 * t267 * t109 * t1rhoa + 0.66725D-1 * t101 * t
     #776 - 0.66725D-1 * t101 * t281 * t783
        t788 = t787 * t294
        t804 = t267 * t103
        t805 = t280 * trhoa
        t815 = cg67 * t100
        t816 = Arhoa * t
        t818 = 0.2D1 * t816 * t1rhoa
        t821 = 0.2D1 * A1rhoa * t * trhoa
        t824 = 0.2D1 * A * t1rhoa * trhoa
        t826 = 0.2D1 * t272 * cg49
        t831 = t275 * t280
        t835 = t280 * t289
        t843 = t101 * t103
        t845 = 0.1D1 / t279 / t107
        t846 = t845 * t289
        t853 = A * t285
        t862 = t104 * t100
        t868 = t815 + t818 + t821 + t824 + t826 + 0.2D1 * A1rhoa * t105
     #* Arhoa + 0.8D1 * t853 * Arhoa * t1rhoa + 0.2D1 * t282 * cg67 + 0.
     #8D1 * t853 * trhoa * A1rhoa + 0.12D2 * t862 * trhoa * t1rhoa + 0.4
     #D1 * t286 * cg49
        t872 = 0.133450D0 * t89 * t1rhoa * t268 + 0.133450D0 * t267 * t7
     #76 * trhoa - 0.133450D0 * t804 * t805 * t783 + 0.133450D0 * t267 *
     # t109 * cg49 + 0.133450D0 * t267 * t276 * t1rhoa + 0.66725D-1 * t1
     #01 * (t815 + t818 + t821 + t824 + t826) * t108 - 0.66725D-1 * t101
     # * t831 * t783 - 0.133450D0 * t804 * t835 * t1rhoa - 0.66725D-1 *
     #t101 * t775 * t280 * t289 + 0.133450D0 * t843 * t846 * t783 - 0.66
     #725D-1 * t101 * t281 * t868
        t875 = t112 ** 2
        t876 = 0.1D1 / t875
        t877 = t293 * t876
        cg47 = -0.2D1 / 0.9D1 * dble(t114) / t298 / t115 * t670
        cg26 = cg78
        t895 = cg7 ** 2
        cg52 = cg44
        t908 = mu ** 2
        t909 = 0.1D1 / t313 / t126 * t908
        t910 = t123 * cg44
        cg65 = 0.2000000001D1 * t315 * cg58 * cg52
        exc_rhoa_rhoa = cg26 * cg25 + cg81 * cg65 + cg78 * cg25 - 0.3D1
     #/ 0.4D1 * rhoa * t7 * cg47 * cg25 + t319 * cg65 + cg81 * cg72 + rh
     #oa * cg26 * cg72 + t129 * (-0.9950248765D1 * t909 * t910 * cg52 +
     #0.2000000001D1 * t315 * cg52 * cg44 + 0.2000000001D1 * t315 * cg58
     # * (cg0 / t304 / cg54 * t121 * t895 + t306 * t310 * cg7 - t306 * t
     #121 * cg47 / 0.2D1 + t120 / t309 / rhoa)) + cg76 + 0.3D1 * t263 *
     #t113 * cg60 + t99 * t788 + cg40 + rho * (cg75 + 0.6D1 * t765 * t26
     #4 * cg60 + 0.3D1 * t263 * t788 * cg18 + 0.3D1 * t263 * t113 * cg30
     # + 0.3D1 * t263 * t295 * cg60 + t99 * t872 * t294 - t99 * t877 * t
     #787)
        cg23 = 0.2D1 * t475
        cg56 = cg45
        t938 = rsrhob * t161
        t945 = t500 * rsrhoa * rsrhob
        t947 = t163 * cg56
        t950 = t164 * rsrhob
        t952 = t13 * cg56
        t954 = rsrhoa * rsrhob
        t956 = t169 * cg56
        cg73 = -0.1328829340D-1 * cg56 * t24 + 0.2137000000D0 * t490 * t
     #334 + 0.2137000000D0 * t938 * t174 - 0.2000000000D1 * t495 * t174
     #* t333 + 0.9999999999D0 * t162 * (-0.1898925000D1 * t945 + 0.37978
     #50000D1 * t947 + 0.35876D1 * cg56 + 0.1228650000D1 * t950 + 0.2457
     #300000D1 * t952 + 0.9858800D0 * t954 + 0.985880D0 * t956) * t173 +
     # 0.1608182432D2 * t520 * t172 * t522 * t333
        t970 = rsrhob * t180
        t995 = rsrhob * t194
        cg22 = (0.4D1 / 0.9D1 * t588 * cg10 * cg24 + 0.4D1 / 0.3D1 * t57
     # * cg23 + 0.4D1 / 0.9D1 * t594 * cg10 * cg24 - 0.4D1 / 0.3D1 * t60
     # * cg23) * dble(t55)
        t1045 = t65 * cg10 * cg24
        t1074 = -0.4D1 * t63 * t214 * cg23 + (-0.638837320D-2 * cg56 * t
     #37 + 0.2054800000D0 * t528 * t344 + 0.2054800000D0 * t970 * t188 -
     # 0.2000000000D1 * t533 * t188 * t343 + 0.1000000000D1 * t181 * (-0
     #.3529725000D1 * t945 + 0.7059450000D1 * t947 + 0.61977D1 * cg56 +
     #0.2524650000D1 * t950 + 0.5049300000D1 * t952 + 0.12503400D1 * t95
     #4 + 0.1250340D1 * t956) * t187 + 0.3216468318D2 * t551 * t186 * t5
     #53 * t343 - cg73) * f * t66 + t218 * frhob * t66 + 0.4D1 * t219 *
     #t373 + t368 * frhoa * t66 + t70 * cg22 * t66 + 0.4D1 * t221 * t373
     # + 0.4D1 * t369 * t223 + 0.4D1 * t371 * t223 + 0.12D2 * t71 * t104
     #5 + 0.4D1 * t71 * t213 * cg23
        cg11 = cg73 + (0.375735750D-1 * cg56 * t50 - 0.1112500000D0 * t5
     #59 * t354 - 0.1112500000D0 * t995 * t202 + 0.2000000000D1 * t564 *
     # t202 * t353 - 0.9999999999D0 * t195 * (-0.2589250000D1 * t945 + 0
     #.5178500000D1 * t947 + 0.36231D1 * cg56 + 0.6601950000D0 * t950 +
     #0.1320390000D1 * t952 + 0.9934200D0 * t954 + 0.993420D0 * t956) *
     #t201 - 0.2960857464D1 * t582 * t200 * t584 * t353) * f * t68 + cg9
     #3 * frhob * t68 - 0.4D1 * t209 * t365 + cg50 * frhoa * t68 + cg62
     #* cg22 * t68 - 0.4D1 * t211 * t365 - 0.4D1 * t361 * t215 - 0.4D1 *
     # t363 * t215 - 0.12D2 * t63 * t64 * t1045 + t1074
        cg68 = -t657 * cg10 * cg24 / 0.9D1 + t226 * cg23 / 0.3D1 - t662
     #* cg10 * cg24 / 0.9D1 - t228 * cg23 / 0.3D1
        cg15 = cg87
        t1109 = t238 * t246 * cg38 / 0.2D1
        t1112 = t85 * t705 * cg80 / 0.2D1
        cg61 = t682 * t683 * cg38 + t686 * t683 * cg80 / 0.2D1 + t692 -
     #t238 * t87 * cg68 / 0.2D1 + t686 * t696 * cg38 / 0.2D1 + t702 * t6
     #96 * cg80 + t708 - t85 * t243 * (-t674 * cg3 * t76 * cg36 / 0.2D1
     #+ t235 * cg15 * t7) / 0.2D1 + t1109 + t1112 + t719
        cg33 = 0.133450D0 * t727 * t729 * t390 - 0.66725D-1 * t251 * (-c
     #g11 * t89 * t93 + 0.3D1 * t252 * t387 + 0.3D1 * t385 * t256 - 0.12
     #D2 * t90 * t747 * cg38 + 0.3D1 * t90 * t255 * cg68) * t95 - 0.6672
     #5D-1 * t251 * t259 * t390 * t95
        t1163 = cg33 * t100
        t1165 = 0.2D1 * t816 * trhob
        t1166 = Arhob * t
        t1168 = 0.2D1 * t1166 * trhoa
        t1171 = 0.2D1 * A * trhob * trhoa
        t1173 = 0.2D1 * t272 * cg61
        t1184 = t403 * t280
        t1207 = t1163 + t1165 + t1168 + t1171 + t1173 + 0.2D1 * Arhob *
     #t105 * Arhoa + 0.8D1 * t853 * Arhoa * trhob + 0.2D1 * t282 * cg33
     #+ 0.8D1 * t853 * trhoa * Arhob + 0.12D2 * t862 * trhoa * trhob + 0
     #.4D1 * t286 * cg61
        t1211 = 0.133450D0 * t89 * trhob * t268 + 0.133450D0 * t267 * t4
     #04 * trhoa - 0.133450D0 * t804 * t805 * t411 + 0.133450D0 * t267 *
     # t109 * cg61 + 0.133450D0 * t267 * t276 * trhob + 0.66725D-1 * t10
     #1 * (t1163 + t1165 + t1168 + t1171 + t1173) * t108 - 0.66725D-1 *
     #t101 * t831 * t411 - 0.133450D0 * t804 * t835 * trhob - 0.66725D-1
     # * t101 * t1184 * t289 + 0.133450D0 * t843 * t846 * t411 - 0.66725
     #D-1 * t101 * t281 * t1207
        exc_rhoa_rhob = cg28 + cg40 + rho * (cg11 + 0.6D1 * t765 * t264
     #* cg38 + 0.3D1 * t263 * t416 * cg18 + 0.3D1 * t263 * t113 * cg68 +
     # 0.3D1 * t263 * t295 * cg38 + t99 * t1211 * t294 - t99 * t877 * t4
     #15)
        cg16 = 0.2D1 * t150 + 0.2D1 * t475
        cg57 = cg56
        t1222 = t333 ** 2
        t1226 = rsrhob ** 2
        t1227 = t500 * t1226
        t1229 = t163 * cg57
        t1232 = t163 * t1226
        t1234 = t13 * cg57
        t1237 = t169 * cg57
        cg66 = -0.1328829340D-1 * cg57 * t24 + 0.4274000000D0 * t938 * t
     #334 - 0.2000000000D1 * t495 * t1222 * t173 + 0.9999999999D0 * t162
     # * (-0.1898925000D1 * t1227 + 0.3797850000D1 * t1229 + 0.35876D1 *
     # cg57 + 0.1228650000D1 * t1232 + 0.2457300000D1 * t1234 + 0.985880
     #0D0 * t1226 + 0.985880D0 * t1237) * t173 + 0.1608182432D2 * t520 *
     # t1222 * t522
        cg12 = cg51
        t1250 = t343 ** 2
        t1272 = t353 ** 2
        cg85 = cg50
        t1290 = cg24 ** 2
        cg70 = (0.4D1 / 0.9D1 * t588 * t1290 + 0.4D1 / 0.3D1 * t57 * cg1
     #6 + 0.4D1 / 0.9D1 * t594 * t1290 - 0.4D1 / 0.3D1 * t60 * cg16) * d
     #ble(t55)
        f1rhob = frhob
        t1312 = cg85 * f
        t1315 = cg62 * f1rhob
        t1332 = cg83 - cg12
        t1339 = t1332 * f
        t1342 = t70 * f1rhob
        t1351 = -0.4D1 * t63 * t214 * cg16 + (-0.638837320D-2 * cg57 * t
     #37 + 0.4109600000D0 * t970 * t344 - 0.2000000000D1 * t533 * t1250
     #* t187 + 0.1000000000D1 * t181 * (-0.3529725000D1 * t1227 + 0.7059
     #450000D1 * t1229 + 0.61977D1 * cg57 + 0.2524650000D1 * t1232 + 0.5
     #049300000D1 * t1234 + 0.12503400D1 * t1226 + 0.1250340D1 * t1237)
     #* t187 + 0.3216468318D2 * t551 * t1250 * t553 - cg66) * f * t66 +
     #t368 * f1rhob * t66 + 0.4D1 * t369 * t373 + t1332 * frhob * t66 +
     #t70 * cg70 * t66 + 0.4D1 * t371 * t373 + 0.4D1 * t1339 * t373 + 0.
     #4D1 * t1342 * t373 + 0.12D2 * t71 * t65 * t1290 + 0.4D1 * t71 * t2
     #13 * cg16
        cg79 = cg66 + (0.375735750D-1 * cg57 * t50 - 0.2225000000D0 * t9
     #95 * t354 + 0.2000000000D1 * t564 * t1272 * t201 - 0.9999999999D0
     #* t195 * (-0.2589250000D1 * t1227 + 0.5178500000D1 * t1229 + 0.362
     #31D1 * cg57 + 0.6601950000D0 * t1232 + 0.1320390000D1 * t1234 + 0.
     #9934200D0 * t1226 + 0.993420D0 * t1237) * t201 - 0.2960857464D1 *
     #t582 * t1272 * t584) * f * t68 + cg50 * f1rhob * t68 - 0.4D1 * t36
     #1 * t365 + cg85 * frhob * t68 + cg62 * cg70 * t68 - 0.4D1 * t363 *
     # t365 - 0.4D1 * t1312 * t365 - 0.4D1 * t1315 * t365 - 0.12D2 * t63
     # * t618 * t1290 + t1351
        cg13 = cg12 + t1312 * t68 + t1315 * t68 - t367 + t1339 * t66 + t
     #1342 * t66 + t375
        cg71 = -t657 * t1290 / 0.9D1 + t226 * cg16 / 0.3D1 - t662 * t129
     #0 / 0.9D1 - t228 * cg16 / 0.3D1
        cg17 = cg38
        t1364 = cg36 ** 2
        cg43 = cg80
        t1370 = t2 * cg38
        t1379 = t2 * cg80
        cg31 = t682 * t1370 * cg17 + t686 * t1370 * cg43 / 0.2D1 + t1109
     # - t238 * t87 * cg71 / 0.2D1 + t686 * t1379 * cg17 / 0.2D1 + t702
     #* t1379 * cg43 + t1112 - t85 * t243 * (-t674 * t1364 * t76 / 0.2D1
     # + t235 * cg15 * t7) / 0.2D1 + t238 * t246 * cg17 / 0.2D1 + t85 *
     #t705 * cg43 / 0.2D1 + t719
        t1rhob = -t238 * t87 * cg17 / 0.2D1 - t85 * t243 * cg43 / 0.2D1
     #- t247 / 0.2D1
        t1400 = cg13 * t89
        t1402 = t255 * cg17
        t1405 = -t1400 * t93 + 0.3D1 * t90 * t1402
        cg41 = 0.133450D0 * t727 * t390 * t728 * t1405 - 0.66725D-1 * t2
     #51 * (-cg79 * t89 * t93 + 0.3D1 * t385 * t1402 + 0.3D1 * t1400 * t
     #387 - 0.12D2 * t90 * t746 * cg38 * cg17 + 0.3D1 * t90 * t255 * cg7
     #1) * t95 - 0.66725D-1 * t251 * t390 * t1405 * t95
        A1rhob = -0.66725D-1 * t251 * t1405 * t95
        t1439 = A1rhob * t100
        t1441 = 0.2D1 * t272 * t1rhob
        t1442 = t1439 + t1441
        t1443 = t1442 * t108
        t1450 = t1439 + t1441 + 0.2D1 * t282 * A1rhob + 0.4D1 * t286 * t
     #1rhob
        t1454 = 0.133450D0 * t267 * t109 * t1rhob + 0.66725D-1 * t101 *
     #t1443 - 0.66725D-1 * t101 * t281 * t1450
        t1455 = t1454 * t294
        t1471 = t280 * trhob
        t1481 = cg41 * t100
        t1483 = 0.2D1 * t1166 * t1rhob
        t1486 = 0.2D1 * A1rhob * t * trhob
        t1489 = 0.2D1 * A * t1rhob * trhob
        t1491 = 0.2D1 * t272 * cg31
        t1499 = t280 * t411
        t1507 = t845 * t411
        t1527 = t1481 + t1483 + t1486 + t1489 + t1491 + 0.2D1 * A1rhob *
     # t105 * Arhob + 0.8D1 * t853 * Arhob * t1rhob + 0.2D1 * t282 * cg4
     #1 + 0.8D1 * t853 * trhob * A1rhob + 0.12D2 * t862 * trhob * t1rhob
     # + 0.4D1 * t286 * cg31
        t1531 = 0.133450D0 * t89 * t1rhob * t397 + 0.133450D0 * t267 * t
     #1443 * trhob - 0.133450D0 * t804 * t1471 * t1450 + 0.133450D0 * t2
     #67 * t109 * cg31 + 0.133450D0 * t267 * t404 * t1rhob + 0.66725D-1
     #* t101 * (t1481 + t1483 + t1486 + t1489 + t1491) * t108 - 0.66725D
     #-1 * t101 * t1184 * t1450 - 0.133450D0 * t804 * t1499 * t1rhob - 0
     #.66725D-1 * t101 * t1442 * t280 * t411 + 0.133450D0 * t843 * t1507
     # * t1450 - 0.66725D-1 * t101 * t281 * t1527
        t1534 = t415 * t876
        cg89 = -0.2D1 / 0.9D1 * dble(t114) / t419 / t131 * t670
        cg90 = cg69
        t1552 = cg2 ** 2
        cg88 = cg8
        t1565 = 0.1D1 / t434 / t142 * t908
        t1566 = t139 * cg8
        cg92 = 0.2000000001D1 * t436 * cg84 * cg88
        exc_rhob_rhob = cg90 * cg20 + cg32 * cg92 + cg69 * cg20 - 0.3D1
     #/ 0.4D1 * rhob * t7 * cg89 * cg20 + t440 * cg92 + cg32 * cg42 + rh
     #ob * cg90 * cg42 + t145 * (-0.9950248765D1 * t1565 * t1566 * cg88
     #+ 0.2000000001D1 * t436 * cg88 * cg8 + 0.2000000001D1 * t436 * cg8
     #4 * (cg1 / t425 / cg19 * t137 * t1552 + t427 * t431 * cg2 - t427 *
     # t137 * cg89 / 0.2D1 + t136 / t430 / rhob)) + cg13 + 0.3D1 * t263
     #* t113 * cg17 + t99 * t1455 + cg28 + rho * (cg79 + 0.6D1 * t765 *
     #t394 * cg17 + 0.3D1 * t263 * t1455 * cg38 + 0.3D1 * t263 * t113 *
     #cg71 + 0.3D1 * t263 * t416 * cg17 + t99 * t1531 * t294 - t99 * t15
     #34 * t1454)
        t1590 = t237 * t86
        t1592 = t84 * t242
        t1594 = t445 * t150
        cg91 = -t1590 * t683 / 0.2D1 - t1592 * t696 / 0.2D1 - t1594 / 0.
     #2D1
        t1599 = t89 * cg86
        t1602 = t101 * A
        t1603 = cg86 * t108
        t1616 = t816 * cg86
        t1617 = t451 * trhoa
        t1618 = t272 * cg91
        t1627 = t835 * cg86
        t1630 = t450 * A
        t1651 = 0.133450D0 * t1599 * t268 + 0.266900D0 * t1602 * t1603 *
     # trhoa - 0.133450D0 * t804 * t805 * t459 + 0.133450D0 * t267 * t10
     #9 * cg91 + 0.133450D0 * t267 * t276 * cg86 + 0.66725D-1 * t101 * (
     #0.2D1 * t1616 + 0.2D1 * t1617 + 0.2D1 * t1618) * t108 - 0.66725D-1
     # * t101 * t831 * t459 - 0.133450D0 * t804 * t1627 - 0.133450D0 * t
     #1630 * t1627 + 0.133450D0 * t843 * t846 * t459 - 0.66725D-1 * t101
     # * t281 * (0.2D1 * t1616 + 0.2D1 * t1617 + 0.2D1 * t1618 + 0.8D1 *
     # t853 * Arhoa * cg86 + 0.12D2 * t862 * trhoa * cg86 + 0.4D1 * t286
     # * cg91)
        exc_rhoa_norm_drho = cg48 + rho * (0.3D1 * t263 * t464 * cg18 +
     #t99 * t1651 * t294 - t99 * t877 * t463)
        cg64 = -t1590 * t1370 / 0.2D1 - t1592 * t1379 / 0.2D1 - t1594 /
     #0.2D1
        t1677 = t1166 * cg86
        t1678 = t451 * trhob
        t1679 = t272 * cg64
        t1688 = t1499 * cg86
        t1711 = 0.133450D0 * t1599 * t397 + 0.266900D0 * t1602 * t1603 *
     # trhob - 0.133450D0 * t804 * t1471 * t459 + 0.133450D0 * t267 * t1
     #09 * cg64 + 0.133450D0 * t267 * t404 * cg86 + 0.66725D-1 * t101 *
     #(0.2D1 * t1677 + 0.2D1 * t1678 + 0.2D1 * t1679) * t108 - 0.66725D-
     #1 * t101 * t1184 * t459 - 0.133450D0 * t804 * t1688 - 0.133450D0 *
     # t1630 * t1688 + 0.133450D0 * t843 * t1507 * t459 - 0.66725D-1 * t
     #101 * t281 * (0.2D1 * t1677 + 0.2D1 * t1678 + 0.2D1 * t1679 + 0.8D
     #1 * t853 * Arhob * cg86 + 0.12D2 * t862 * trhob * cg86 + 0.4D1 * t
     #286 * cg64)
        exc_rhob_norm_drho = cg48 + rho * (0.3D1 * t263 * t464 * cg38 +
     #t99 * t1711 * t294 - t99 * t1534 * t463)
        t1717 = cg86 ** 2
        t1721 = A * t1717
        t1726 = t280 * cg86 * t459
        t1732 = t459 ** 2
        t1746 = t463 ** 2
        exc_norm_drho_norm_drho = rho * (t99 * (0.133450D0 * t89 * t1717
     # * t109 + 0.667250D0 * t101 * t1721 * t108 - 0.266900D0 * t804 * t
     #1726 - 0.266900D0 * t1630 * t1726 + 0.133450D0 * t101 * t103 * t84
     #5 * t1732 - 0.66725D-1 * t101 * t281 * (0.2D1 * t1721 + 0.12D2 * t
     #862 * t1717)) * t294 - t99 * t1746 * t876)
        exc_rhoa_norm_drhoa = cg81 * cg5 + t319 * cg5 + t129 * (-0.99502
     #48765D1 * t909 * t910 * cg4 + 0.2000000001D1 * t315 * cg4 * cg44 +
     # 0.2000000001D1 * t315 * cg58 * (-t305 * t121 * cg7 / 0.2D1 - t119
     # * t310 / 0.2D1))
        t1766 = cg4 ** 2
        exc_norm_drhoa_norm_drhoa = t129 * (-0.9950248765D1 * t909 * t12
     #3 * t1766 + 0.2000000001D1 * t315 * t1766)
        exc_rhob_norm_drhob = cg32 * cg39 + t440 * cg39 + t145 * (-0.995
     #0248765D1 * t1565 * t1566 * cg74 + 0.2000000001D1 * t436 * cg74 *
     #cg8 + 0.2000000001D1 * t436 * cg84 * (-t426 * t137 * cg2 / 0.2D1 -
     # t135 * t431 / 0.2D1))
        t1790 = cg74 ** 2
        exc_norm_drhob_norm_drhob = t145 * (-0.9950248765D1 * t1565 * t1
     #39 * t1790 + 0.2000000001D1 * t436 * t1790)
        cg94 = exc_norm_drhob_norm_drhob
        return
      end"""

def maple2f90(code,replacements={}):
    c1=code.replace("\n     #","")

    floatRe=re.compile(r"([0-9]+\.?[0-9]*|\.[0-9]+)[dD]([-+]?[0-9]+)")
    c1=floatRe.sub(r"\1e\2_dp",c1)

    varNameRe=re.compile(r"([a-zA-Z][a-zA-Z0-9_]*)")
    c2=varNameRe.split(c1)
    for i in range(len(c2)):
        if c2[i] in replacements.keys():
            c2[i]=replacements[c2[i]]
    c1="".join(c2)
    vars=re.findall(r" *([a-zA-Z][a-zA-Z0-9_]*) *=",c1)
    lc1=["real(kind=dp) :: "+", ".join(vars)+"\n","\n"]+c1.splitlines()
    splittableRe=re.compile(r"([ /=()])")
    for i in range(len(lc1)):
        if (len(lc1[i])>80):
            pieces=splittableRe.split(lc1[i])
            ll=0
            newL=[]
            for p in pieces:
                if ll<10 or ll+len(p)<75:
                    newL.append(p)
                    ll+=len(p)
                else:
                    newL.append("&\n      ")
                    newL.append(p)
                    ll=len(p)+6
            lc1[i]="".join(newL)
    c1="\n".join(lc1)+"\n"
    return c1

#===============================================================================
if __name__=="__main__":
    replacements={}
    for i in range(len(renamedVar)):
        replacements[renamedVar[i]]=origNames[i]
    output = maple2f90(code,replacements)

    # print output unless we are selftesting
    print(output)
