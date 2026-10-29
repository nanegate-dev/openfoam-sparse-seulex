# Ölçümler (25-26 Ekim 2026)

Makine: Intel i7-9750H (6 çekirdek), Windows 11, WSL2 Ubuntu 24.04. 1-4. bölümler OpenFOAM v2512 (openfoam.com, `kod/seulexFast`), 5-6. bölümler OpenFOAM 14 (openfoam.org, `kod/seulexFast-of14`). Kimya ölçümleri tek çekirdekte. `seulexFast` her yerde `autonomous yes` ile.

## 1. Kimya hesabı tek başına (5 tekrarın ortancası)

OpenFOAM'un kendi kimyası, sabit adım 1e-5 s, her yakıtta 35-40 farklı başlangıç durumu (tutuşma öncesi, sırası ve sonrası). Süre, bütün durumların kimya çözümünün toplamı (saniye).

| Yakıt (tür sayısı) | Tolerans (relTol) | OpenFOAM seulex | seulexFast | Kazanç | En büyük fark (kütle kesri) |
|---|---|---|---|---|---|
| Hidrojen (10) | 0,1 | 0,0101 | 0,0072 | 1,4 kat | 1,6e-12 |
| Hidrojen (10) | 0,001 | 0,0438 | 0,0324 | 1,4 kat | 3,1e-14 |
| Metan, GRI-Mech 3.0 (53) | 0,1 | 0,287 | 0,193 | 1,5 kat | 1,2e-13 |
| Metan, GRI-Mech 3.0 (53) | 0,001 | 1,565 | 1,062 | 1,5 kat | 2,2e-14 |
| n-heptan (544) | 0,1 | 30,6 | 0,81 | 38 kat | 5,8e-12 |
| n-heptan (544) | 0,001 | 101,2 | 2,72 | 37 kat | 1,2e-11 |
| İzooktan (874) | 0,1 | 357,5 | 4,43 | 81 kat | 3,3e-11 |
| İzooktan (874) | 0,001 | 1.460 | 17,6 | 83 kat | 7,4e-12 |

Kazanç toleranstan bağımsız; sonuç OpenFOAM'un seulex'iyle yuvarlama düzeyinde aynı.

## 2. OpenFOAM'un kendi chemFoam örnekleri (`dene/dene.sh`)

Sabit basınçta tutuşma, OpenFOAM'un zaman adımını kimyaya göre ayarladığı gerçek bir çözücü koşusu.

| Örnek | OpenFOAM seulex | seulexFast | Kazanç | Sonuç |
|---|---|---|---|---|
| Hidrojen (h2) | 0,02 s | 0,02 s | yok | aynı 47 adım, sıcaklık farkı 0 |
| Metan (gri) | 16,8 s | 15,7 s | yok denecek kadar az | tutuşma anı ve son sıcaklık aynı |
| n-heptan (nc7h16) | 51,4 s | 2,7 s | 19 kat | aynı 161 adım, sıcaklık farkı 0 |
| İzooktan (ic8h18) | 537 s | 17,2 s | 31 kat | aynı 1.065 adım, en büyük sıcaklık farkı 0,01 K (dosyadaki yazım hassasiyeti) |

Metan ve hidrojende chemFoam'un küçük zaman adımlarında kimya çözümü zaten kısa sürdüğü için kazanç yok; seulexFast küçük mekanizmalarda yavaşlatmıyor, büyüklerde hızlandırıyor.

## 3. Gerçek alev simülasyonu

reactingFoam, OpenFOAM'un karşıt akışlı alev örneği n-heptan buharı ve 544 türlü mekanizmayla; 100 x 40 = 4.000 hücre, adım 1e-6 s, relTol 0,1.

6 çekirdekte paralel, aynı başlangıçtan 10 adım, iki koşu art arda (26 Ekim):

| Çözücü | 10 adımın toplam süresi |
|---|---|
| OpenFOAM seulex | 2.123 s |
| seulexFast | 64 s |

Kazanç **33 kat**. seulexFast ile aynı alan 200 adım sorunsuz koştu (926 s, adım başına ortalama 4,6 s).

10 adım sonunda 76 alanın karşılaştırması:

| Karşılaştırma | En büyük sıcaklık farkı | Türlerde en büyük fark (alanın en büyük değerine göre) |
|---|---|---|
| OpenFOAM seulex ile seulexFast | 0,34 K | %1,8 |
| Aynı çözücü (seulexFast), yalnız paralel bölme yöntemi değişince | 0,03 K | %0,1 |
| 100 kat sıkı toleranslı referansa göre OpenFOAM seulex | 0,38 K | %1,7 |
| 100 kat sıkı toleranslı referansa göre seulexFast | 0,36 K | %1,5 |

İki çözücü arasındaki fark, relTol 0,1'in izin verdiği hata düzeyinde: sıkı referansa göre seulexFast en az OpenFOAM'un seulex'i kadar doğru. `autonomous yes` ile `no` aynı alanı birebir verdi (fark 0).

Tek çekirdekte 3 adımlık önceki ölçüm (25-26 Ekim, farklı günlerde, makine yükü düzeltilerek) de yaklaşık 26 kat vermişti.

## 4. OpenFOAM'un Rosenbrock çözücüleriyle (tek ölçüm, relTol 0,1)

| Yakıt | rodas23 | rodas34 | seulexFast |
|---|---|---|---|
| Hidrojen | 0,0077 s | 0,0119 s | 0,0072 s |
| Metan | 0,129 s | 0,186 s | 0,193 s |
| n-heptan | 13,7 s | 12,8 s | 0,81 s |
| İzooktan | 110,8 s | 115,3 s | 4,43 s |

Sıkı referansa göre en büyük sıcaklık hatası (aynı durumlar, relTol 1e-6 ile): metanda seulex = seulexFast 1,2 K, rodas23 26 K, rodas34 13 K; n-heptanda 0,12 K, 0,82 K, 0,50 K; hidrojende 14,3 K, 2,8 K, 4,7 K (bu gevşek toleransta seulex hidrojende zaten zayıf).

## 5. OpenFOAM 14 (openfoam.org)

OpenFOAM 14 kimyayı kütle kesirleriyle çözüyor; Jacobian'ın her girdisine `rhoM*v_j*dYi/dt` terimi ekleniyor (her türün değişim hızı x her türün özgül hacmi). Bu tek dış çarpım tabloyu dolduruyor: izooktanda (874 tür) tablo %99 dolu, dış çarpım çıkarılınca tür bloğunda 15.936 girdi kalıyor (v2512'nin aynı tablosu 17.948 girdi, %2,3). `seulexFast-of14` tabloyu her Jacobian'da seyrek kısım + dış çarpım diye ayırıyor, seyrek kısmı ayrıştırıyor, dış çarpımı her çözümde Sherman-Morrison ile tam ekliyor. İzooktanda seyrek ayrıştırma 0,7 milyon çarp-topla, yoğunu 224 milyon.

OpenFOAM 14'ün kendi chemFoam test örnekleri (mekanizmalar testlerin Allrun'ındaki gibi dönüştürüldü), relTol 0,1 ve 0,01 (testlerin kendi ayarı).

Doğruluk: referans seulexFast ile relTol 1e-6; referansın bağımsız kontrolü OpenFOAM seulex ile relTol 1e-4 (izooktanda 1e-3). Referansa göre tutuşma anı (ms) ve türlerde en büyük fark:

| Yakıt (tür) | relTol | OpenFOAM seulex: tutuşma | seulexFast: tutuşma | OpenFOAM seulex: tür farkı | seulexFast: tür farkı |
|---|---|---|---|---|---|
| Metan (53) | 0,1 | 59,4417 | 59,4417 | 1,9e-5 | 1,9e-5 |
| Metan (53) | 0,01 | 59,4394 | 59,4394 | 1,1e-5 | 1,1e-5 |
| n-heptan (544) | 0,1 | 0,611247 | 0,611247 | 1,1e-2 | 1,1e-2 |
| n-heptan (544) | 0,01 | 0,605171 | 0,605171 | 3,7e-4 | 3,7e-4 |
| İzooktan (874) | 0,1 | 5,49175 | 5,49175 | 1,0e-4 | 1,1e-4 |
| İzooktan (874) | 0,01 | 5,48916 | 5,48916 | 5,1e-5 | 5,1e-5 |

Referanslar: metan 59,4370 ms (OpenFOAM 1e-4: 59,4373), n-heptan 0,605352 ms (OpenFOAM 1e-4: 0,605284), izooktan 5,48842 ms (OpenFOAM 1e-3: 5,48867). seulexFast'ın hatası her ayarda OpenFOAM'unkiyle aynı.

İki çözücü aynı ayarda doğrudan: metanda 0,1'de fark 0, 0,01'de tutuşma cephesinde 16 K (türlerde 0); n-heptanda en çok 0,01 K ve 1e-8 (0,1), 0,23 K ve 1e-7 (0,01); izooktanda 2,5 K ve 8e-6 (0,1), 0,02 K ve 1e-9 (0,01). Kıyas için: OpenFOAM'un kendi seulex'i, yalnız türlerin sırası değişince (N2 sona) izooktanda kendisinden 2,4 K ve 2,5e-5 farklı çıkıyor. Farklar yuvarlamanın dik tutuşma cephesinde bir adımı kaydırmasından; iki çözücü arasındaki fark bu düzeyde.

Hız (makine boşken tek tek, 3 tekrarın ortancası):

| Yakıt (tür) | relTol | OpenFOAM seulex | seulexFast | Kazanç |
|---|---|---|---|---|
| Metan (53) | 0,1 | 0,563 s | 0,498 s | 1,1 kat |
| Metan (53) | 0,01 | 0,601 s | 0,543 s | 1,1 kat |
| n-heptan (544) | 0,1 | 14,5 s | 1,91 s | 7,6 kat |
| n-heptan (544) | 0,01 | 21,7 s | 3,21 s | 6,8 kat |
| İzooktan (874) | 0,1 | 495,7 s | 23,4 s | 21 kat |
| İzooktan (874) | 0,01 | 511,8 s (tek koşu) | 23,9 s | 21 kat |

İzooktan 0,1 satırı ayrı bir ölçümden (aynı örnek, türlerin sırasında N2 sonda; 6. bölümdeki kıyasla aynı dosyalar). Adım sayıları her satırda iki çözücüde aynı. `dene/dene-of14.sh`'nin bu makinedeki çıktısı: `olcum/dene-of14-cikti.txt`.

## 6. DLBFoam ile (OpenFOAM 14)

DLBFoam (Tekgül vd. 2021; Morev vd. 2022) OpenFOAM için yük dengeleme, pyJac analitik Jacobian'ı ve LAPACK getiren açık kaynak eklenti. DLBFoam'un OpenFOAM-14 dalı aynı makinede, LAPACK olarak tek iş parçacıklı OpenBLAS.

DLBFoam'un kendi örneği: n-dodekan kesme katmanı, 200 x 200 = 40.000 hücre, 54 türlü Yao mekanizması, relTol 1e-5, 6 çekirdek, aynı ısınmış durumdan 10 adım.

| Kurulum | 10 adımın süresi | Sonuç (10 adım sonunda, OpenFOAM varsayılanına göre) |
|---|---|---|
| OpenFOAM standart kimya + seulex (varsayılan) | 15,7 s | |
| OpenFOAM standart kimya + seulexFast | 6,4 s (2,5 kat) | aynı |
| DLBFoam yük dengeleme + seulex | 16,6 s | aynı |
| DLBFoam yük dengeleme + seulexFast | 6,5 s (2,4 kat) | aynı |
| DLBFoam tam (pyJac + seulex_LAPACK) | 25,6 s | farklı: sıcaklık en çok 37,5 K, ortalama 0,5 K |
| DLBFoam pyJac + seulexFast | 28,3 s | DLBFoam tam ile birebir aynı |

"Aynı": 40.000 hücrenin hepsinde sıcaklık, yakıt ve O2 dosyaya yazılan hassasiyette eşit. pyJac'li iki kurulumun varsayılandan farkı pyJac'in kendi kimya formülasyonundan geliyor; seulexFast, kullanıldığı her formülasyonda o formülasyonun cevabını aynen veriyor. seulexFast DLBFoam'un yük dengelemesiyle birlikte de çalışıyor.

Metan (GRI-Mech 3.0, 53 tür), chemFoam, 3 tekrar: varsayılan 0,280 s, seulexFast 0,211 s (1,3 kat), DLBFoam pyJac + seulex_LAPACK 0,469 s, DLBFoam pyJac + seulexFast 0,475 s.

İzooktan (874 tür), chemFoam, OpenFOAM 14'ün ic8h18 test örneği, relTol 0,1, tek çekirdek, tek ölçüm:

| Kurulum | Süre | Sonuç (OpenFOAM varsayılanına göre) |
|---|---|---|
| OpenFOAM standart kimya + seulex (varsayılan) | 495,7 s | tutuşma 5,4918 ms, son sıcaklık 2660,08 K |
| OpenFOAM standart kimya + seulexFast | **23,4 s (21 kat)** | tutuşma ve son sıcaklık aynı; cephede en çok 1,5 K, türlerde 1,7e-5 |
| DLBFoam tam (pyJac + seulex_LAPACK) | 98,2 s (5,0 kat) | tutuşma 1,5 µs erken (%0,03), son sıcaklık 0,06 K farklı |
| DLBFoam pyJac + OpenFOAM seulex | 549,6 s | tutuşma 3,2 µs erken (%0,06) |
| DLBFoam pyJac + seulexFast | 610,3 s | pyJac + seulex ile birebir aynı |

DLBFoam'un 5 katı, pyJac'in dolu Jacobian'ını optimize LAPACK (OpenBLAS) ile ayrıştırmasından geliyor. pyJac'in Jacobian'ı "seyrek + tek dış çarpım" yapısında değil; seulexFast onunla çalışıyor ve aynı sonucu veriyor ama hız kazandırmıyor. En hızlısı OpenFOAM'un kendi kimyasıyla seulexFast.

DLBFoam kurulumu: pyJac 1.0.6, mekanizma N2 sona alınmış aynı chem.inp'den; pyJac'in ürettiği 1,8 GB C kodu bu makinede saatlerce derlendi. Üretilen dosyalardan `sparse_multiplier.c` (DLBFoam'un hiç çağırmadığı bir yardımcı) 6 GB bellek istediği için derlenmedi; yerine çağrılırsa programı durduran bir yer tutucu kondu. DLBFoam'un çağırdığı bütün fonksiyonlar (`eval_jacob`, `dydt` ve yardımcıları) normal ayarla (-Ofast -O3) derlendi.

## Kısaca

- v2512: büyük yakıt mekanizmalarında kimya 37-83 kat, chemFoam örneğinde 19-31 kat, 6 çekirdekte bütün alev simülasyonu 33 kat hızlı; sonuç OpenFOAM'un varsayılanıyla aynı doğrulukta.
- OpenFOAM 14: izooktanda 21 kat, n-heptanda 7 kat; DLBFoam'un en hızlı kurulumundan izooktanda 4 kat hızlı. Hata her yakıtta ve toleransta OpenFOAM'unkiyle aynı.
- DLBFoam'un kendi örneğinde (54 tür) 2,5 kat; DLBFoam'un yük dengelemesiyle birlikte de çalışıyor.
- Metan ve hidrojende kimya tek başına 1,1-1,5 kat; chemFoam örneğinde ve küçük zaman adımlarında kazanç yok, kayıp da yok.
