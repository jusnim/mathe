# Klausurstoff Diskrete Mathematik: Was ist am wichtigsten?

Kurs: Diskrete Mathematik, HTW Berlin.
Stand: 28.09.2026.
Quellen: nur die Dateien unter `latex/wissen/` und `latex/uebungen/`.

Dieses Dokument hat fünf Teile:

1. Methode: Wie wurde gewichtet?
2. Themen, sortiert nach Wahrscheinlichkeit
3. Übungsaufgaben, sortiert nach Priorität
4. Plan für Phase 3: Python-Skripte
5. Fehler und Unsicherheiten in den Quellen

Kurze Erklärung der Abkürzungen:

| Abkürzung | Bedeutung | Datei |
|---|---|---|
| K | echte Klausur vom 29.07.2026 (Foto) | `uebungen/92_klausur_pzr1.tex` |
| GP | Gedächtnisprotokoll (Klausur oder Probeklausur) | `uebungen/90_gedaechtnisprotokoll.tex` |
| KV | Klausurvorbereitung 1–8 mit Lösungen | `uebungen/91_klausurvorbereitung.tex` |
| Ü1 … Ü6 | Übungsblatt 1 … 6 mit Lösungen | `uebungen/0X_uebung.tex` |
| VL1 … VL6 | Vorlesung 1 … 6 | `wissen/0X_*.tex` |
| Z. | Zeile in der Datei | – |

Hinweis zu den Fundstellen: Die Sätze und Definitionen in den Vorlesungsdateien haben **keine Nummern**.
Deshalb steht bei jeder Fundstelle: Datei, Abschnitt, Name der Umgebung und Zeilennummer.

---

## 1. Methode: Wie wurde gewichtet?

### 1.1 Grundidee

Eine Aufgabe, die schon in einer Klausur stand, kommt sehr wahrscheinlich wieder.
Eine Aufgabe, die nur in der Vorlesung steht, kommt selten.
Deshalb zählt jeder Beleg unterschiedlich viel.

### 1.2 Gewicht der Belege

| Beleg | Gewicht pro Vorkommen | Begründung |
|---|---|---|
| K (echte Klausur 29.07.2026) | 5 | Sicherer Beleg. So sieht die Klausur aus. |
| GP (Gedächtnisprotokoll) | 3 | Wahrscheinlich eine frühere Klausur. Die Aufgaben sind aber nur sinngemäß erinnert. |
| KV (Klausurvorbereitung mit Lösungen) | 2 | Der Dozent hat diese Aufgaben selbst für die Klausur ausgewählt. |
| Ü1–Ü6 (Übungsblätter) | 1 | Zeigen die Aufgabenformen. Viele Themen, aber nicht alle sind klausurrelevant. |
| VL1–VL6 (Vorlesung) | 0,5 | Nur Theorie. Zeigt, wo das Wissen steht. |

### 1.3 Von Punkten zu Prozent

Für jedes Thema wurden die Gewichte addiert.
Danach wurde die Summe in eine **geschätzte** Wahrscheinlichkeit übersetzt.
Diese Prozentzahlen sind keine Statistik.
Sie sind eine begründete Schätzung.

Grobe Regel:

| Belege | Geschätzte Wahrscheinlichkeit |
|---|---|
| in K **und** GP **und** KV | 90–95 % |
| in K und KV, aber nicht in GP | 75–85 % |
| nur als Teilschritt einer Klausuraufgabe | 70–90 % |
| nur in KV und Übung | 25–40 % |
| nur in Übung | 5–20 % |
| nur in der Vorlesung | unter 5 % |

### 1.4 Wichtigste Beobachtung

K und GP haben **fast dieselbe Struktur**:

| Nr. in K | Nr. in GP | Thema |
|---|---|---|
| 1 | 1 | RSA (verschlüsseln, entschlüsseln, Oscar berechnet d) |
| 2 | 2 | Inklusion-Exklusion (Zahlen ≤ N, nicht teilbar durch …) |
| 3 | 3 | Euklidischer Algorithmus und diophantische Gleichung |
| 4 | – | Chinesischer Restsatz (baut auf Aufgabe 3 auf) |
| 5 | 4 | [10, 8, 3]₁₁-Code: Kontrollmatrix, Syndrom, Decodieren |

Die empfangenen Wörter in GP Aufgabe 4 sind fast gleich wie in K Aufgabe 5.
Die RSA-Zahlen in GP Aufgabe 1 sind dieselben wie in Ü4 Aufgabe 9.
Folgerung: Der Dozent nutzt ein **festes Klausurschema**.
Die Zahlen ändern sich, die Aufgabentypen bleiben.

Bewertung in K: 29 Punkte = 1,0; 15 Punkte = 4,0 (Datei K, Z. 77–80).

---

## 2. Themen, sortiert nach Wahrscheinlichkeit

Begriffe kurz erklärt:

- **ggT**: größter gemeinsamer Teiler.
- **Euklidischer Algorithmus (EA)**: wiederholte Division mit Rest. Der letzte Rest ungleich 0 ist der ggT.
- **Lemma von Bezout**: Es gibt ganze Zahlen x, y mit a·x + b·y = ggT(a, b).
- **Matrix Q**: Produkt der Matrizen (qᵢ 1; 1 0) aus dem EA. Ihre Determinante liefert x und y.
- **Diophantische Gleichung**: Gleichung, deren Lösungen ganze Zahlen sein müssen.
- **φ(m) (Eulersche Phi-Funktion)**: Anzahl der Zahlen von 1 bis m, die mit m teilerfremd sind.
- **Syndrom**: S(r) = H·rᵀ. Es zeigt, ob ein empfangenes Wort r ein Codewort ist.
- **Kontrollmatrix H**: Matrix mit H·cᵀ = 0 genau für Codewörter c.
- **Erzeugermatrix G**: Matrix, deren Zeilen den Code erzeugen. Codieren: c = a·G.

| Rang | Thema | Wahrsch. | Belege (Datei + Aufgabe) | Fundstelle im Wissen | Typische Aufgabenform |
|---|---|---|---|---|---|
| 1 | **RSA-Verfahren** | 95 % | K A1; GP A1; KV A7; Ü4 A9 | VL3 `wissen/03_potenzen_zm.tex` Abschn. 3.3: Definition „RSA: Ver- und Entschlüsselung“ (Z. 160), Beispiel m = 101 (Z. 172), Oscar-Beispiel (Z. 196), Definition „Öffentlicher und privater Schlüssel“ (Z. 207); `wissen/91_klausurvorbereitung.tex` Bemerkung „RSA: Berechnung des privaten Schlüssels d“ (Z. 32) | e binär schreiben; c = wᵉ mod m berechnen; Formel w = cᵈ mod m nennen; φ(m) berechnen; d auf zwei Wegen berechnen (EA/Bezout und Euler) |
| 2 | **Euklidischer Algorithmus, ggT, lineare diophantische Gleichung** | 95 % | K A3; GP A3; KV A3; Ü2 A4 | VL1 `wissen/01_ring_ganze_zahlen.tex`: Lemma von Bezout (Z. 406), Folgerung Lösbarkeit (Z. 411), Division mit Rest (Z. 439), EA (Z. 445–471), Satz Euklid (Z. 473), Beweis mit Matrix Q (Z. 481), Beispiel 1965/225 (Z. 494); `wissen/91_klausurvorbereitung.tex` Algorithmus EA (Z. 12) | a) ggT(a, b) mit EA; b) eine Lösung von a·x − b·y = c; c) alle Lösungen |
| 3 | **[10, 8, 3]₁₁-Code (modifizierter ISBN-Code, „Modulo-11-Code“)** | 90 % | K A5; GP A4; Ü6 6_9; Ü5 5_8 (Parameter) | VL6 `wissen/06_decodierung.tex` Abschn. 6.3 (Z. 215–281): H, systematische Form, Formeln für c₉ und c₁₀ (Z. 246), Syndrom und Fehlerstelle x = s₂/s₁ (Z. 253–260), Bemerkung „s₁ = 0, s₂ ≠ 0“ (Z. 279); `wissen/06_uebung.tex` Inverse mod 11 (Z. 28), Algorithmus Decodierung (Z. 38) | H aufschreiben; Kontrollziffern c₉, c₁₀ berechnen; Syndrome berechnen; entscheiden, welches Wort nicht decodierbar ist; das andere decodieren; Probe H·cᵀ = 0 |
| 4 | **Inklusion-Exklusion (Siebformel)** | 90 % | K A2; GP A2; KV A6; Ü4 A1 | VL2 `wissen/02_zaehlprinzipien.tex` Abschn. 2.5: Formel r = 3 mit α₁, α₂, α₃ (Z. 189–198), Satz „Allgemein“ (Z. 200), Beispiel Musikschule (Z. 210); `wissen/91_klausurvorbereitung.tex` Siebformel (Z. 25) | Wie viele Zahlen ≤ N sind durch keine der Zahlen p₁, p₂, p₃ teilbar? Mit Aᵢ und αᵢ aufschreiben. Oder: Vereinsmitglieder mit Sportarten |
| 5 | **Schnelles Potenzieren (Square-and-Multiply)** | 90 % (als Teil von RSA) | K A1a,b; KV A7 (1. Weg); Ü4 A9; Ü4 A7 | VL3 `wissen/03_potenzen_zm.tex` Abschn. 3.1: Algorithmus (Z. 15), Beispiel 12¹⁰⁰ mod 34 (Z. 22), Bemerkung k mod φ(m) (Z. 129), Beispiel 11¹⁶ mod 34 (Z. 137) | Exponent binär schreiben; Tabelle a, a², a⁴, … mod m; passende Potenzen multiplizieren |
| 6 | **Primfaktorzerlegung und φ(n)** | 85 % | K A1d; GP A1 („Primzahlzerlegung“); KV A4; Ü3 A2; Ü4 A9 (φ(φ(m))) | VL1 `wissen/01_ring_ganze_zahlen.tex` Definition φ (Z. 639), Satz „Berechnung von φ(m)“ (Z. 643); VL2 `wissen/02_zaehlprinzipien.tex` Satz „Formel für φ(n)“ (Z. 226), Beispiel φ(60) (Z. 256); VL1 Hauptsatz der Arithmetik (Z. 419) | n in Primfaktoren zerlegen; φ(n) = n·∏(1 − 1/p) oder φ(n) = ∏φ(pᵏ) |
| 7 | **Chinesischer Restsatz (simultane Kongruenzen)** | 80 % | K A4; KV A8; Ü3 A1 | VL1 `wissen/01_ring_ganze_zahlen.tex` Satz Chinesischer Restsatz (Z. 580), Beweis „Bestimmung von f⁻¹“ (Z. 589), Beispiel mit 3 Kongruenzen (Z. 603); `wissen/91_klausurvorbereitung.tex` Simultane Kongruenzen (Z. 44) | z ≡ a mod m, z ≡ b mod n; über m·x − n·y = b − a lösen; alle Lösungen z₀ + k·m·n; kleinste positive Lösung |
| 8 | **Inverse in ℤₘ (über Bezout oder Euler)** | 75 % (als Teilschritt) | K A1e, K A5e (Division mod 11); Ü2 A5; Ü4 A8; Ü6 6_9 (Inversentabelle) | VL1 Satz „[a] invertierbar ⇔ ggT(a, m) = 1“ (Z. 298); `wissen/02_uebung.tex` Bemerkung Bezout/Inverse (Z. 17); `wissen/04_uebung.tex` Satz a⁻¹ = a^(φ(m)−1) (Z. 36); `wissen/06_uebung.tex` Tabelle Inverse mod 11 (Z. 28) | a⁻¹ mod m mit EA rückwärts; Tabelle der Inversen mod 11 |
| 9 | **ISBN-10: Prüfziffer und Zahlendreher** | 35 % | KV A1, KV A2; Ü2 A1, Ü2 A2 | VL1 `wissen/01_ring_ganze_zahlen.tex` Abschn. „Der ISBN Code“ (Z. 341–400): Prüfsumme Σ i·cᵢ ≡ 0 (Z. 376–381), Satz „ISBN10 erkennt genau einen Fehler“ (Z. 388); `wissen/02_uebung.tex` Bemerkung Prüfziffer (Z. 10) | c₁₀ ergänzen (X = 10); zeigen, dass ein Zahlendreher erkannt wird |
| 10 | **Lineare Codes allgemein: G, H, systematische Form, Syndromdecodierung, Hamming-Codes** | 30 % | Ü5 5_1 bis 5_5; Ü6 6_4, 6_8; (Grundlage für K A5) | VL4 `wissen/04_lineare_codes.tex`: Erzeugermatrix (Z. 137), Kontrollmatrix (Z. 174), systematische Form (Z. 181–184), Dualitätssatz (Z. 201), binärer Hamming-Code (Z. 217), Kodieren (Z. 248), Dekodier-Algorithmus (Z. 275), q-närer Hamming-Code (Z. 293), Klassenführer (Z. 336, Algorithmus Z. 341); VL6 Syndromdecodierung (`wissen/06_decodierung.tex` Z. 122) | G und H angeben; Wort codieren; prüfen, ob Codewort; Syndrom → Fehlerstelle → korrigieren |
| 11 | **Turnierplan (Rundenturnier)** | 25 % | KV A5; Ü1 A3, Ü1 A4 | VL1 `wissen/01_ring_ganze_zahlen.tex` Algorithmus Rundenturnier (Z. 157), Beispiele 2m = 4 und 6 (Z. 146, 163); `wissen/91_klausurvorbereitung.tex` Bemerkung Turnierplan (Z. 21) | Plan für 2m Mannschaften aufstellen; Gegner von Mannschaft 2m in Runde r |
| 12 | **Schranken: Hamming, Singleton; perfekt, optimal (MDS)** | 20 % | Ü5 5_6, 5_7, 5_8 | VL4 `wissen/04_lineare_codes.tex` Abschn. 4.5: Hamming-Schranke (Z. 386), perfekt (Z. 390), Singleton-Schranke (Z. 416), MDS-Code (Z. 423); `wissen/05_uebung.tex` (Z. 27, 33) | Schranke für [n, k, d]_q prüfen; ist der Code perfekt oder optimal? kleinstes n bestimmen |
| 13 | **Fermat, Euler, Ordnung, primitive Elemente** | 20 % | Ü4 A7, Ü4 A8; KV A7 (Satz von Euler als Weg) | VL3 `wissen/03_potenzen_zm.tex` Satz Fermat (Z. 70), Definition Ordnung/primitiv (Z. 80), Satz Euler (Z. 110) | aᵏ mod p mit Fermat; Ordnung von Elementen in ℤₚ*; primitive Elemente |
| 14 | **Einheitengruppe ℤₘ*, Nullteiler, Gruppentafel** | 10 % | Ü2 A3 | VL1 Einheiten (Z. 273–339); VL3 Bemerkung ℤₘ* (Z. 43) | Einheiten und Inverse auflisten; Gruppentafel |
| 15 | **Standardarray, Nebenklassen, Klassenführer** | 15 % | Ü6 6_8 | VL6 `wissen/06_decodierung.tex` Abschn. 6.1: Nebenklasse (Z. 25), Satz (Z. 31), Klassenführer/Standardarray (Z. 58), Algorithmus (Z. 76) | Klassenführer und Syndrome tabellieren; Wörter decodieren |
| 16 | **Zyklische Codes: Generatorpolynom, Kreisteilungsklassen, X^n − 1 zerlegen** | 15 % | Ü6 6_1, 6_2, 6_3, 6_5, 6_6, 6_7 | VL5 `wissen/05_zyklische_codes.tex`: Definition (Z. 9), Satz Generatorpolynom (Z. 162), Generatormatrix (Z. 217), Kodieren (Z. 250), Syndrom (Z. 261), Kreisteilungsklassen (Z. 284), Kontrollpolynom (Z. 353), reziprokes Polynom (Z. 421); VL6 Syndrom als Polynomrest (`wissen/06_decodierung.tex` Z. 137–213); `wissen/06_uebung.tex` (Z. 10, 14) | X^n − 1 in irreduzible Polynome zerlegen; g(X) und zyklische G angeben |
| 17 | **Derangements (fixpunktfreie Permutationen)** | 10 % | Ü4 A2 | VL2 `wissen/02_zaehlprinzipien.tex` Beispiel Derangement (Z. 267) | d₄ berechnen und auflisten |
| 18 | **Teilbarkeitsregeln (9, 11, 7)** | 10 % | Ü1 A1, Ü1 A2 | VL1 Satz Teilbarkeitsregeln (Z. 119) | Quersumme, alternierende Quersumme, Regel mit 1001 |
| 19 | **Kombinatorik: Binomialkoeffizient, Gitterwege, binomischer Lehrsatz** | 10 % | Ü3 A4–A7 | VL2 Kombinationen (Z. 44), Kombinationen mit Wiederholung (Z. 75), Binomischer Lehrsatz (Z. 145), Folgerungen (Z. 156) | Wege im Gitter zählen; Summen mit Binomialkoeffizienten |
| 20 | **Teileranzahl d(n), Teilersumme σ(n), vollkommene Zahlen, Mersenne-Primzahlen** | 5 % | Ü4 A3–A6 | `wissen/04_uebung.tex` Definition (Z. 8), Satz Formeln (Z. 18), vollkommen (Z. 26), Mersenne (Z. 30) | Formeln herleiten; 6, 28, 496 prüfen |
| 21 | **BIBD (Blockpläne)** | 5 % | Ü1 A5 | VL1 Definition Parameter BIBD (Z. 177), Satz (Z. 189) | b·k = v·r und r(k−1) = λ(v−1) beweisen |
| 22 | **Affine Ebene AG(2, p)** | 5 % | Ü3 A3 | keine Fundstelle in `wissen/` | Punkte einer Geraden; Schnittpunkt |
| 23 | Gruppenaxiome, Homomorphismen, Grinberg-Formel, Partitionen | unter 5 % | keine Übungsaufgabe | VL1 Gruppe (Z. 9), Homomorphismus (Z. 524), Grinberg (Z. 208); VL2 Partitionen (Z. 309) | nur Theorie |

---

## 3. Übungsaufgaben, sortiert nach Priorität

Die Klausuraufgaben (K, GP) stehen oben.
Sie haben keine Lösung.
Sie sind die wichtigsten Testfälle.

„Lösung“ heißt: Die Datei enthält eine Lösung des Dozenten.

| Rang | Datei + Aufgabe | Thema | Lösung | Kurzbeschreibung |
|---|---|---|---|---|
| 1 | K A1 (`92_klausur_pzr1.tex`) | RSA | nein | m = 43·57 = 2451, e = 221, w = 1511: e binär, c, Entschlüsselungsformel, φ(m), d |
| 2 | K A5 | [10,8,3]₁₁-Code | nein | a = 11500005: H, c₉, c₁₀; r₁ = 1150000511, r₂ = 1150000733 decodieren |
| 3 | K A3 | EA, Diophantik | nein | ggT(17, 15); Lösung von 17x − 15y = 1 |
| 4 | K A2 | Inklusion-Exklusion | nein | Zahlen ≤ 720, weder durch 3, 4 noch 5 teilbar |
| 5 | K A4 | Chinesischer Restsatz | nein | z ≡ 8 mod 17, z ≡ 5 mod 15; mit Ergebnis aus A3; kleinste positive Lösung |
| 6 | GP A1 (`90_gedaechtnisprotokoll.tex`) | RSA | nein | m = 2803, „d = 113“, w = 715; Formeln; Oscar: zwei Wege |
| 7 | GP A4 | [10,8,3]₁₁-Code | nein | H; Syndrome von r₁, r₂; nicht decodierbares Wort; korrigieren; Probe |
| 8 | GP A2 | Inklusion-Exklusion | nein | Zahlen ≤ 840, nicht durch 2, 3, 5 teilbar |
| 9 | GP A3 | EA, Diophantik | nein | a·x − b·y = 30; Werte unbekannt |
| 10 | KV A7 (`91_klausurvorbereitung.tex`) | RSA | ja | d für m = 101, e = 13 auf zwei Wegen (Euler, EA/Q) |
| 11 | Ü4 A9 (`04_uebung.tex`) | RSA | ja | m = 2803, e = 113, a = 715: c = 708; d = 1463 auf zwei Wegen; Probe |
| 12 | Ü6 6_9 (`06_uebung.tex`) | [10,8,3]₁₁-Code | ja | r₁, r₂, r₃ decodieren mit x = s₁/s₀ |
| 13 | KV A3 | EA, Diophantik | ja | ggT(2406, 654); 2406x − 654y = 24; alle Lösungen |
| 14 | Ü2 A4 (`02_uebung.tex`) | EA, Diophantik | ja | 1001x + 840y = 98 und = 0; alle Lösungen (Matrix Q) |
| 15 | KV A6 | Inklusion-Exklusion | ja | 55 Athleten, 3 Sportarten: wie viele machen etwas anderes? |
| 16 | Ü4 A1 | Inklusion-Exklusion | ja | Zahlen 1 … 1000, nicht durch 2, 3, 5 teilbar (266) |
| 17 | KV A8 | Chinesischer Restsatz | ja | z ≡ 5 mod 13, z ≡ 4 mod 15 (109) |
| 18 | Ü3 A1 (`03_uebung.tex`) | Chinesischer Restsatz | ja | zwei Systeme; EA rückwärts und Matrix Q |
| 19 | KV A4 | Primfaktoren, φ | ja | 360 zerlegen; φ(360) = 96 auf zwei Wegen |
| 20 | Ü3 A2 | φ | ja | φ(1000), φ(1001), φ(1002) |
| 21 | Ü2 A5 | Inverse in ℤₘ | ja | 6⁻¹ in ℤ₁₁ und ℤ₁₇, 3⁻¹ in ℤ₁₀, 5⁻¹ in ℤ₁₂ |
| 22 | Ü4 A7 | Schnelles Potenzieren, Fermat | ja | 3¹⁰⁰⁰ mod 7 |
| 23 | KV A1 | ISBN-10 | ja | Prüfziffer von 3-05-501517-? |
| 24 | Ü2 A1 | ISBN-10 | ja | drei Prüfziffern (7, X, 7) |
| 25 | KV A2 | ISBN-10 Zahlendreher | ja | zeigen: h·rᵀ ≠ 0 |
| 26 | Ü2 A2 | ISBN-10 Zahlendreher | ja | wie KV A2 |
| 27 | Ü5 5_8 (`05_uebung.tex`) | Schranken | ja | [10,8,3]₁₁ und [10,6,5]₁₁: Hamming- und Singleton-Schranke |
| 28 | Ü4 A8 | Ordnung, primitiv, Inverse | ja | Ordnung von 2, 3, 4, 5 in ℤ₁₇* |
| 29 | Ü6 6_8 | Standardarray | ja | H₃(2): Klassenführer, Syndrome, r₁ = 2200, r₂ = 0121 decodieren |
| 30 | Ü5 5_4 | Hamming H₂(3) | ja | drei Wörter decodieren |
| 31 | Ü5 5_5 | Hamming H_q(m) | ja (mit Fehler) | Parameter; H für H₃(3); r₁, r₂ decodieren |
| 32 | Ü5 5_1 | Wiederholungscode | ja | G, H; codieren; Codewörter prüfen |
| 33 | Ü5 5_2 | Paritätscode | ja | G, H; codieren; Codewörter prüfen |
| 34 | Ü5 5_3 | Hamming H₂(2), H₂(4) | ja | H und G angeben |
| 35 | Ü6 6_4 | Hamming H₃(2) | ja | Parameter, G, H |
| 36 | KV A5 | Turnierplan | ja | 8 Mannschaften, 7 Spieltage |
| 37 | Ü1 A3 (`01_uebung.tex`) | Turnierplan | ja | wie KV A5 |
| 38 | Ü1 A4 | Turnierplan | ja | Gegner von Mannschaft 2m in Runde r |
| 39 | Ü5 5_6 | Hamming-Schranke | ja | Golay-Codes [23,12,7]₂, [11,6,5]₃ sind perfekt |
| 40 | Ü5 5_7 | Schranken | ja (unvollständig) | kleinstes n für [n,5,5]₂ |
| 41 | Ü2 A3 | Einheitengruppe | ja | ℤ₁₂*, ℤ₁₈*, Inverse, Gruppentafeln |
| 42 | Ü6 6_1 | Zyklische Codes | ja | irreduzible Polynome Grad 3; X⁷ − 1 zerlegen |
| 43 | Ü6 6_2 | Zyklische Codes | ja (mit Fehler) | g(X) und zyklische G für H₂(3) |
| 44 | Ü6 6_3 | Zyklische Codes | ja | g(X) und G für [7,6,2]₂ |
| 45 | Ü6 6_5 | Zyklische Codes | ja | X⁴ − 1 in F₃[X] zerlegen |
| 46 | Ü6 6_6 | Zyklische Codes | ja | zwei ternäre zyklische [4,2]-Codes |
| 47 | Ü6 6_7 | Zyklische Codes | ja | Gibt es einen zyklischen Code isomorph zu H₃(2)? |
| 48 | Ü4 A2 | Derangements | ja | fixpunktfreie Permutationen in S₄ |
| 49 | Ü1 A1 | Teilbarkeit | ja | 299792 durch 9 oder 11 teilbar? |
| 50 | Ü1 A2 | Teilbarkeit | ja | Regel für 7 mit 1001 |
| 51 | Ü3 A4 | Kombinatorik | ja | Gitterwege im 4×4-Gitter (70) |
| 52 | Ü3 A5 | Kombinatorik | ja | Gitterwege im m×n-Gitter |
| 53 | Ü3 A6 | Binomischer Lehrsatz | ja | alternierende Summe = 1 |
| 54 | Ü3 A7 | Binomischer Lehrsatz | ja | Vandermonde-Identität beweisen |
| 55 | Ü3 A3 | Affine Ebene | ja | Geraden in AG(2,5) |
| 56 | Ü4 A3 | d(n), σ(n) | ja | Formeln herleiten |
| 57 | Ü4 A4 | vollkommene Zahlen | ja | 6 und 28 sind vollkommen |
| 58 | Ü4 A5 | Mersenne | ja | drei Primzahlen 2ᵏ − 1 |
| 59 | Ü4 A6 | Mersenne, vollkommen | ja | 2^(k−1)(2ᵏ − 1) ist vollkommen; 496 |
| 60 | Ü1 A5 | BIBD | ja | Gleichungen für BIBDs beweisen |

---

## 4. Plan für Phase 3: Python-Skripte

### 4.1 Allgemeine Regeln für alle Skripte

- Ein Skript pro **Aufgabentyp**. Nicht pro Aufgabe.
- Die Ausgabe ist ein **Lösungsweg wie auf Papier**. Jeder Schritt steht in einer eigenen Zeile. Die Sprache ist Deutsch.
- Jede Konvention ist eine **Option** (Kommandozeilen-Schalter). Der Standardwert ist die Konvention der Klausur bzw. des Dozenten.
- Jedes Skript hat einen Testmodus. Er rechnet die Übungsaufgaben nach und vergleicht mit den Lösungen des Dozenten.
- **Repräsentanten** (welche Zahl für eine Restklasse steht) sind überall eine Option:
  `--rep standard` (0 … m−1) oder `--rep symmetrisch` (−m/2 … m/2).
  Beleg: Die Quellen nutzen beides. Beispiele: d = −1339 ≡ 1463 (Ü4 A9), d = −23 ≡ 77 (KV A7), z = −86 (KV A8), Inversentabelle mit −5 … 5 (Ü6 6_9), 12⁴ ≡ −4 mod 34 (VL3 Z. 29).
- Die Reihenfolge unten ist: zuerst nach Klausurpriorität, aber Hilfsskripte stehen vor den Skripten, die sie benutzen.

### 4.2 Liste der Skripte

#### 1. `euklid.py` – EA, ggT, Bezout, diophantische Gleichung, Inverse

- **Eingaben:** a, b; optional rechte Seite c; Form der Gleichung (`+` oder `-`); optional Modul m für die Inverse.
- **Ausgabe:**
  1. EA-Tabelle: a = q₀·b + r₁, b = q₁·r₁ + r₂, … bis Rest 0.
  2. ggT = letzter Rest ≠ 0.
  3. Lösbarkeit: Teilt ggT die Zahl c?
  4. Eine Lösung, auf zwei Wegen: (a) „EA rückwärts“ mit allen Einsetzschritten; (b) Matrix Q = Q₀·…·Qₙ mit Zwischenprodukten und Det Q = ±1.
  5. Hochmultiplizieren mit c/ggT.
  6. Homogene Lösung und allgemeine Lösung mit Parameter t.
  7. Probe durch Einsetzen.
  8. Modus Inverse: a⁻¹ mod m mit Probe a·a⁻¹ ≡ 1.
- **Testfälle:** K A3 (17, 15; 17x − 15y = 1); KV A3 (2406x − 654y = 24; Lösung x = 3, y = 11 und x = 112, y = 412); Ü2 A4 (1001x + 840y = 98; x = 658, y = −784); Ü2 A5 (Inverse 6 mod 11 = 2, 6 mod 17 = 3, 3 mod 10 = 7, 5 mod 12 = 5); VL1 Beispiel 1965/225 (Z. 449, 494); GP A3 (Form a·x − b·y = 30, Werte selbst wählen).
- **Konventionen, die zu verschiedenen Lösungen führen:**

| Punkt | Varianten | Beleg in den Quellen |
|---|---|---|
| Vorzeichenform der Gleichung | a·x + b·y = c oder a·x − b·y = c | „+“: Ü2 A4, VL1 Lemma (Z. 406). „−“: K A3, KV A3, GP A3 |
| Lösungsweg | EA rückwärts oder Matrix Q | beide in KV A3 b) und Ü3 A1; nur Q in Ü2 A4 und Ü4 A9 |
| Vorzeichen von Det Q | Det Q = (−1)^(n+1); Faktor ±r_n | VL1 Z. 487–490; Beispiel 1965/225 mit Det Q = −1 und Faktor −15 (Z. 501) |
| Welche spezielle Lösung | direkt aus dem EA (z. B. 17·(−7) − 15·(−8) = 1) oder die mit kleinstem positiven x (x = 8, y = 9) | KV A3: EA rückwärts gibt (3, 11), Matrix Q gibt (112, 412) – beide richtig |
| Richtung des Parameters t | x = x₀ + (b/g)·t, y = y₀ + (a/g)·t (bei „−“) oder x = x₀ − (b/g)·t, y = y₀ + (a/g)·t (bei „+“) | KV A3 c): x = 3 + 109t, y = 11 + 401t; Ü2 A4 (iii): (658, −784) + (−120, 143)·t |
| Reihenfolge der Eingaben im EA | EA mit (a, b) oder mit (b, a) | Ü3 A1 a): „EA für 7, 5“ bei z = 5x + 3 = 7y + 4; KV A7: „EA mit φ(m) und e“ |
| Repräsentant der Inversen | 0 … m−1 oder symmetrisch | Ü2 A5 c): 3⁻¹ = −3 = 7; Ü6 6_9 Tabelle mit −5 … 5 |

#### 2. `primfaktor_phi.py` – Primfaktorzerlegung und φ(n)

- **Eingaben:** n; optional Liste bekannter Faktoren.
- **Ausgabe:** Probedivision Schritt für Schritt; Zerlegung n = ∏ pᵏ; φ(n) auf zwei Wegen: (a) n·∏(1 − 1/p), (b) ∏(pᵏ − p^(k−1)); optional φ(φ(n)).
- **Testfälle:** KV A4 (360 → 96); Ü3 A2 (1000 → 400, 1001 → 720, 1002 → 332); VL2 Beispiel φ(60) = 16 (Z. 256); Ü4 A9 (φ(2802) = 932); K A1 d (2451 = 3·19·43 → φ = 1512).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Rechenweg | Produktformel n·∏(1 − 1/p) oder Multiplikativität | KV A4 zeigt beide („oder so“); Ü3 A2 nur Multiplikativität |
| Angabe „m = p·q“ | Faktoren als Primzahlen annehmen oder **immer** vollständig zerlegen | K A1: m = 43·57, aber 57 = 3·19 ist **keine** Primzahl. Das Skript muss immer vollständig zerlegen und warnen |

#### 3. `schnell_potenzieren.py` – aᵏ mod m, Fermat/Euler, Ordnung

- **Eingaben:** a, k, m; Schalter für Reduktion mit Euler/Fermat; Modus „Ordnung“ (alle Potenzen bis 1).
- **Ausgabe:** k binär; Tabelle a^(2^i) mod m; Auswahl der Potenzen mit Bit 1; Produkt Schritt für Schritt mod m; optional k mod φ(m) vorab (nur wenn ggT(a, m) = 1); Ordnung und „primitiv ja/nein“.
- **Testfälle:** VL3 12¹⁰⁰ mod 34 = 30 (Z. 22); VL3 11¹⁶ mod 34 = 1 (Z. 137); VL3 8¹³ mod 101 = 18 (Z. 172); Ü4 A7 (3¹⁰⁰⁰ mod 7 = 4); Ü4 A8 (Ordnungen in ℤ₁₇*: 8, 16, 4, 16); Ü4 A9 Tabellen (715^(2^i) mod 2803, 113^(2^i) mod 2802).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Schreibweise der Binärzahl | höchstes Bit links (1110001 für 113) oder als Summe von Zweierpotenzen | Ü4 A9: „binär e = 1110001“; VL3: 100 = 2⁶ + 2⁵ + 2²; KV A7: 39 = 1 + 2 + 4 + 32 |
| Repräsentanten in der Tabelle | 0 … m−1 oder symmetrisch | VL3 Beispiel 12¹⁰⁰: −4, −16; Ü4 A8: 2, 4, 8, −1, … |
| Reduktion des Exponenten | nur binär oder zuerst k mod φ(m) | VL3 Bemerkung Z. 129 (nur bei ggT(a, m) = 1; Warnbeispiel 12 mod 34, Z. 121) |

#### 4. `rsa.py` – RSA komplett

- **Eingaben:** m (oder Faktoren); e **oder** d; Klartext w bzw. a; oder Codewort c.
- **Ausgabe:**
  1. e binär (nutzt `schnell_potenzieren.py`).
  2. Formel c = wᵉ mod m und Rechnung.
  3. Formel w = cᵈ mod m.
  4. Primfaktorzerlegung von m und φ(m) (nutzt `primfaktor_phi.py`).
  5. Weg 1 für d: e·d ≡ 1 mod φ(m) mit EA und Matrix Q (nutzt `euklid.py`).
  6. Weg 2 für d: d ≡ e^(φ(φ(m)) − 1) mod φ(m) mit schnellem Potenzieren.
  7. Probe: cᵈ mod m = w.
  8. Text „Was braucht Oscar?“: m faktorisieren → φ(m) → d.
- **Testfälle:** K A1 (m = 2451, e = 221, w = 1511 → e = 11011101₂, c = 896, φ = 1512, d = 821); Ü4 A9 (m = 2803, e = 113, a = 715 → c = 708, d = 1463); KV A7 (m = 101, e = 13 → d = 77); VL3 Beispiel (m = 101, e = 13, a = 8 → c = 18); VL3 Bob-Beispiel (p = 101, q = 113, e = 3533 → d = 6597, a = 9726 → c = 5761); GP A1 (m = 2803, „d = 113“, w = 715).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| m prim oder zusammengesetzt | φ(m) = m − 1 oder Zerlegung | Ü4 A9 und KV A7: m prim. K A1: m = 43·57 zusammengesetzt. GP A1 Notiz: „zwei cases, m Primzahl und m nicht Primzahl“ |
| φ(m) bei m = p·q | (p−1)(q−1) nur, wenn p und q prim sind | K A1: falsch wäre (43−1)(57−1) = 2352 → d = 149 (entschlüsselt falsch). Richtig: φ = 42·2·18 = 1512 → d = 821 |
| gegebener Schlüssel | e gegeben oder d gegeben | GP A1 nennt „d = 113“; Ü4 A9 mit denselben Zahlen nennt e = 113. Falls d = 113 wirklich privat ist: e = 1463 und c = 715¹⁴⁶³ mod 2803 = 265 statt 708. Beide Varianten ausgeben |
| Name des Klartexts | a oder w | K A1: w; Ü4 A9 und VL3: a |
| Weg für d | EA/Matrix Q oder Euler-Potenz | KV A7 und Ü4 A9 zeigen beide; GP A1: „zwei Lösungswege angeben“ |
| Repräsentant von d | negativ aus Det Q oder positiv | Ü4 A9: d = −1339 ≡ 1463; KV A7: d = −23 ≡ 77 |

#### 5. `mod11_code.py` – [10, 8, 3]₁₁-Code (und Zusatz ISBN-Varianten)

- **Eingaben:** Klartextwort a (8 Ziffern) **oder** empfangenes Wort r (10 Ziffern); Optionen unten.
- **Ausgabe:**
  1. H in Grundform und (optional) in systematischer Form (A | E).
  2. Formeln c₉ ≡ Σ(1+i)·cᵢ und c₁₀ ≡ Σ(9−i)·cᵢ, jede Summe ausgeschrieben; Codewort c.
  3. Syndrom S(r) = H·rᵀ mit ausgeschriebenen Summen.
  4. Fallunterscheidung: S = 0 → Codewort. s₁ ≠ 0 und x = s₂/s₁ ∈ {1, …, 10} → ein Fehler. Sonst → nicht decodierbar (mindestens 2 Fehler), mit Begründung.
  5. Fehlerstelle x (Division mod 11 über Inversentabelle), Fehlergröße, korrigiertes c.
  6. Probe H·cᵀ = 0.
  7. Warnung, wenn r nicht 10 Stellen hat (GP A4).
- **Testfälle:** K A5 (a = 11500005 → c₉ = 4, c₁₀ = 6, c = 1150000546; r₁ → S = (3, 0) → nicht decodierbar; r₂ → S = (9, 10) → x = 6, Fehlergröße 9, c = 1150020733); Ü6 6_9 (r₁ → x = 3, y = 4; r₂ → x = 6, y = −2; r₃ Codewort); VL6 Beispiele (r = 0…050 → x = 9, e = 5; r = 35…000 → x = 3, e = 8; `wissen/06_decodierung.tex` Z. 262, 270); GP A4 (Wörter mit fehlender Stelle).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Form von H | Grundform (1 1 … 1 ; 1 2 … 10) oder systematisch (A \| E) | VL6 Z. 221–229: beide. Ü6 6_9: Grundform. Achtung: Die Syndrome sind verschieden. Für K r₁: Grundform (3, 0), systematisch (8, 6). Für r₂: (9, 10) bzw. (3, 6). Standard: Grundform |
| Reihenfolge der Zeilen von H | Einsenzeile oben oder Zeile 1 … 10 oben | VL6 und Ü6 6_9: Einsenzeile oben. Option für umgekehrte Reihenfolge |
| Lage der Einheitsmatrix | H = (A \| E) mit E rechts (VL6 6.3) oder H = (E \| A) mit E links (VL4, Ü6 6_4) | VL6 Z. 224: „Einheitsmatrix rechts“; VL4 Z. 182: H = (E \| A) |
| Namen der Syndromteile | (s₁, s₂) oder (s₀, s₁) | VL6 Z. 257: s₁, s₂. Ü6 6_9 und `wissen/06_uebung.tex` Z. 51: s₀, s₁ |
| Syndrom als Spalte oder Zeile | Spaltenvektor (H·rᵀ) oder Zeile | alle Quellen: Spalte. Ü6 6_8 schreibt in der Tabelle als Zeile „1 0“ |
| Nummerierung der Stellen | 1 … 10 oder 0 … 9 | VL1, VL6, Ü6 6_9, K A5 (c₉, c₁₀): ab 1. Mit Zählung ab 0 wird x falsch |
| Repräsentanten | 0 … 10 oder −5 … 5 | Ü6 6_9: H mit −5 … −1, Inverse mit −5 … 5 |
| Kriterium „nicht decodierbar“ | nur „s₁ = 0 und s₂ ≠ 0“ (VL6 Z. 279) **oder zusätzlich** „s₂ = 0 und s₁ ≠ 0“ (x = 0 ist keine Stelle) | K r₁ hat genau den zweiten Fall. Die Vorlesung nennt ihn nicht ausdrücklich |
| Kodierung | c = a·G mit G = (E \| −Aᵀ) oder direkte Formeln für c₉, c₁₀ | VL6 Z. 233–248: beide sind gleich; Ausgabe beider Wege als Kontrolle |

#### 6. `inklusion_exklusion.py` – Siebformel

- **Eingaben:** Modus „Teiler“: N und Liste von Teilern. Modus „Mengen“: |Ω|, |Aᵢ|, |Aᵢ ∩ Aⱼ|, … direkt. Modus „Derangement“: n.
- **Ausgabe:** Definition Ω und Aᵢ; jede Schnittmenge |A_i ∩ A_j| = ⌊N / kgV⌋; α₁, α₂, α₃; Formel |Ω| − α₁ + α₂ − α₃; Ergebnis; Kontrolle mit Produktformel und Hinweis, ob sie exakt ist.
- **Testfälle:** K A2 (N = 720; 3, 4, 5 → 720 − 564 + 144 − 12 = 288); GP A2 (N = 840; 2, 3, 5 → 224); Ü4 A1 (N = 1000; 2, 3, 5 → 266); KV A6 (55 Athleten → 4); VL2 Musikschule (73 Kinder → 11); Ü4 A2 (d₄ = 9).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Schnittmengen | ⌊N/(a·b)⌋ oder ⌊N/kgV(a, b)⌋ | K A2 enthält 4 (keine Primzahl). Hier sind 3, 4, 5 paarweise teilerfremd, also gleich. Bei Teilern wie 2 und 4 nur kgV richtig |
| Abrunden | ⌊N/d⌋ | Ü4 A1: ⌊1000/3⌋ = 333 |
| Produktformel N·∏(1 − 1/p) | nur Kontrolle; exakt nur, wenn N durch das Produkt teilbar ist | Ü4 A1: „nur näherungsweise“ (266,66). K A2: 720 ist durch 60 teilbar → exakt 288 |
| Schreibweise | Aᵢ und αᵢ, Indizes nach Teiler (A₂, A₃, A₅) oder 1, 2, 3 | K A2 und GP A2 fordern „Symbole Aᵢ, αᵢ aus der Vorlesung“; Ü4 A1 nutzt A₂, A₃, A₅; VL2 und KV A6 nutzen A₁, A₂, A₃ |

#### 7. `crt.py` – Chinesischer Restsatz

- **Eingaben:** Liste von Paaren (aᵢ, mᵢ); optional eine vorgegebene Bezout-Lösung (für „mit Hilfe von Aufgabe 3“).
- **Ausgabe:** Ansatz z = m·x + a = n·y + b; Gleichung m·x − n·y = b − a; Lösbarkeit (ggT | b − a); EA (nutzt `euklid.py`); spezielle Lösung; z₀; alle Lösungen z₀ + k·m·n; kleinste positive Lösung; Probe. Bei mehr als zwei Kongruenzen schrittweise zusammenfassen.
- **Testfälle:** K A4 (z ≡ 8 mod 17, z ≡ 5 mod 15 → z ≡ 110 mod 255; über K A3); KV A8 (→ −86 ≡ 109 mod 195); Ü3 A1 a (→ 18 mod 35); Ü3 A1 b (→ −73 ≡ 47 mod 60); VL1 Beispiel mit 3 Kongruenzen (→ −74 ≡ 1246 mod 1320, Z. 603).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Welche Kongruenz ist „m, a“ | erste oder zweite | VL1 Z. 595 und KV A8: erste. Ü3 A1 a: EA „für 7, 5“ |
| Vorzeichen | m·x − n·y = b − a | VL1 Z. 596; KV A8; Ü3 A1 b (b − a = 3) |
| Spezielle Lösung | aus eigener EA-Rechnung oder aus der Vorgabe | K A4 a: „mit Hilfe des Ergebnisses aus Aufgabe 3“. Mit (x, y) = (8, 9) aus 17x − 15y = 1 folgt z = −400; mit (−7, −8) folgt z = 365. Beide ≡ 110 mod 255 |
| Darstellung des Ergebnisses | negativer Repräsentant oder kleinster positiver | KV A8: z = −86, dann +195 = 109; Ü3 A1 b: nur −73 mod 60 |

#### 8. `isbn10.py` – ISBN-10 Prüfziffer und Fehlererkennung

- **Eingaben:** 9 Ziffern (Prüfziffer berechnen) oder 10 Zeichen (prüfen); Modus „Zahlendreher an Stellen x, y“.
- **Ausgabe:** Summe Σ i·cᵢ mit allen Summanden; Reduktion mod 11; X für 10; Prüfung h·cᵀ ≡ 0; beim Zahlendreher h·eᵀ = (x − y)(c_y − c_x).
- **Testfälle:** KV A1 und Ü2 A1 (3-05-501517 → 7; 0-471-82819 → X; 4-263-79215 → 7); VL1 Beispiele 0-387-98999-4, 3-528-06580-X (Z. 364); KV A2, Ü2 A2.
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Gewichte | Σ i·cᵢ mit i = 1 … 10 (Quellen) oder 10 … 1 (Praxis) | VL1 Z. 377–380, Ü2 A1. Beide geben dieselbe Prüfziffer, weil Σ(11 − i)cᵢ ≡ −Σ i·cᵢ mod 11. Die Zwischensummen sind aber verschieden |
| Rechentrick | Gewichte 6 … 9 als −5 … −2 | Ü2 A1 und KV A1: „1.3 + … − 5.1 − 4.5 − 3.1 − 2.7“ |
| Wert 10 | als X schreiben | VL1 Z. 381; Ü2 A1 |

#### 9. `linearer_code.py` – lineare Codes, Hamming-Codes, Syndromdecodierung, Standardarray

- **Eingaben:** q; G oder H (als Matrix) oder Code-Typ (Wiederholung, Parität, H₂(m), H_q(m)); Wort a oder r; Optionen unten.
- **Ausgabe:** Umrechnung G ↔ H in systematischer Form; Parameter [n, k, d]; Codieren a·G; Syndrom H·rᵀ; Vergleich mit Spalten von H (Fehlerstelle x, Fehlergröße y); korrigiertes c; Tabelle Klassenführer/Syndrom; d über den Dualitätssatz.
- **Testfälle:** Ü5 5_1 (Wiederholungscode F₃⁵); Ü5 5_2 (Paritätscode F₃⁵); Ü5 5_3 (H₂(2), H₂(4)); Ü5 5_4 (Decodieren in H₂(3)); Ü5 5_5 (H₃(3); richtige Stelle für r₂ ist x = 11); Ü6 6_4 und 6_8 (H₃(2), Klassenführer, r₁ = 2200 → 2210, r₂ = 0121 → 0221); VL4 Beispiel H₃(2) (r₁ = 1111, r₂ = 0011 → 1011, Z. 351); VL6 Beispiel [4,2,2]₂ Nebenklassen (Z. 44).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Systematische Form | (A) H = (E \| A), G = (−Aᵀ \| E); (B) G = (E \| A), H = (−Aᵀ \| E); (C) H = (A \| E), G = (E \| −Aᵀ) | (A): VL4 Z. 182–183, Ü6 6_4, Ü6 6_8. (B): Ü5 5_1, Ü5 5_2. (C): VL6 Z. 224–233. Alle drei als Option |
| Spaltenreihenfolge von H für H₂(3) | verschiedene Reihenfolgen der Spalten 4–7 | VL4 Z. 229: (110, 011, 111, 101). Ü5 5_4: (110, 101, 011, 111). VL6 Z. 167: wie VL4. Verschiedene H geben verschiedene Syndrome und G. H muss frei wählbar sein |
| Spalten von H_q(m) | Repräsentant mit erster Nicht-Null = 1 („Klassenführer“) oder beliebig | VL4 Z. 337; Ü5 5_5 c (erste Nicht-Null = 1) |
| Kodierung | c = a·G (Information rechts, VL4: c₁ c₂ c₃ a₁ … a₄) oder Information links (Ü5 5_1, 5_2) | VL4 Z. 249; Ü5 5_2 b: a·G = 2 1 2 0 2 |
| Syndrom | Spalte H·rᵀ oder Zeile r·Hᵀ | alle Quellen: H·rᵀ als Spalte. Ausgabe beider Schreibweisen |
| Nummerierung | Stellen ab 1 | VL4 Z. 268 (Stelle x, hₓ = x-te Spalte); Ü5 5_5 (x = 13) |
| Klassenführer bei Gleichstand | welcher Vektor mit kleinstem Gewicht | VL6 Beispiel [4,2,2]₂: 1000, 0100, 0010 gewählt; 0001 wäre auch möglich (`wissen/06_decodierung.tex` Z. 49–54). Option für Auswahlregel |

#### 10. `schranken.py` – Hamming-Schranke, Singleton-Schranke, perfekt, optimal

- **Eingaben:** [n, k, d]_q; oder k, d, q für „kleinstes n“.
- **Ausgabe:** e = ⌊(d − 1)/2⌋; |B(e)| = Σ C(n, i)(q − 1)ⁱ ausgeschrieben; Vergleich |C|·|B(e)| mit qⁿ; perfekt ja/nein; Singleton k + d ≤ n + 1; MDS/optimal ja/nein; Tabelle für n = 9, 10, … bis die Schranke erfüllt ist.
- **Testfälle:** Ü5 5_6 (Golay: perfekt); Ü5 5_7 (Ergebnis fehlt in der Quelle; Rechnung ergibt n ≥ 12 nach Hamming-Schranke); Ü5 5_8 ([10,8,3]₁₁: 101 vs. 121; [10,6,5]₁₁: 4601 vs. 14641); VL4 Beispiel H₂(m) perfekt (Z. 394); VL5 Golay-Beweise.
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Bedeutung von „optimal“ | Singleton mit Gleichheit (MDS) | Ü5 5_8 nennt Codes mit k + d = n + 1 „optimal“; VL4 Z. 423 nennt sie MDS |
| e aus d | 2e + 1 = d bzw. ⌊(d−1)/2⌋ | `wissen/05_uebung.tex` Z. 29 |

#### 11. `turnierplan.py` – Rundenturnier

- **Eingaben:** Anzahl 2m der Mannschaften; optional Runde r.
- **Ausgabe:** Tabelle aller Runden; Regel x + y ≡ r mod (2m − 1); Sonderfall 2x ≡ r → Spiel gegen 2m; Gegner von 2m: x ≡ r·m mod (2m − 1).
- **Testfälle:** KV A5, Ü1 A3 (8 Mannschaften), Ü1 A4 (Gegner von 8), VL1 Beispiele 2m = 4 und 2m = 6 (Z. 146, 163).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Repräsentant der Restklasse 0 | 0 oder 2m − 1 | Ü1 A3: Runden r = 1 … 7; Paar (1, 7) in Runde 1; Ü1 A4: „7 ≡ 7 (mod 7)“ |
| Formel | x + y ≡ r mod (2m − 1) (Ü1, VL1) oder „x + y = 2m“ (KV A5, unklar formuliert) | KV A5 und Ü1 A3 haben dieselbe Tabelle; Standard: Formel aus VL1 Z. 159 |

#### 12. `einheitengruppe.py` – ℤₘ*, Inverse, Gruppentafel, Ordnungen

- **Eingaben:** m; Modus (Einheiten, Inversentabelle, Gruppentafel, Ordnungen).
- **Ausgabe:** Liste der Einheiten mit ggT-Test; Nullteiler; Inversentabelle; Gruppentafel; Ordnung jedes Elements; primitive Elemente.
- **Testfälle:** Ü2 A3 (ℤ₁₂*, ℤ₁₈*); Ü4 A8 (ℤ₁₇*); Ü6 6_9 (Inverse mod 11); VL1 ℤ₆ (Z. 277); VL3 ℤ₇* (Z. 94).
- **Konventionen:** Repräsentanten symmetrisch (Ü2 A3: {±1, ±5}; Ü6 6_9) oder 0 … m−1 (Ü4 A8 Probe: 2·9 ≡ 1).

#### 13. `zyklischer_code.py` – zyklische Codes über F_p

- **Eingaben:** p, n; optional g(X), k, Wort y.
- **Ausgabe:** Kreisteilungsklassen von ℤₙ; Zerlegung von Xⁿ − 1 (irreduzible Faktoren durch Probieren); alle Generatorpolynome; zyklische G (Verschiebungen von g); h(X), reziprokes h, H; Codieren a(X)·g(X); Syndrom als Rest y(X) mod g(X); Decodieren mit Klassenführer-Tabelle.
- **Testfälle:** Ü6 6_1 (X⁷ − 1 über F₂); Ü6 6_2 (H₂(3), zwei g); Ü6 6_3 (g = X + 1); Ü6 6_5 (X⁴ − 1 über F₃); Ü6 6_6, 6_7; VL5 Beispiele n = 4, 5, 9, 15, Golay G₂₃ und G₁₁; VL6 Beispiel y(X) = X⁶ + X + 1 → c(X) = X⁶ + X⁴ + X + 1 (Z. 196–211).
- **Konventionen:**

| Punkt | Varianten | Beleg |
|---|---|---|
| Nummerierung der Stellen | ab 0 (c₀ … c_{n−1}) | VL5 Z. 12, 53. Anders als ISBN und mod-11-Code (ab 1) |
| Reihenfolge Wort ↔ Polynom | a₀ a₁ … ↔ a₀ + a₁X + … (niedrigster Grad links) | VL5 Z. 53; Ü6 6_2 (g = 1 + X + X³ → 1101000); VL6 Z. 208–211 |
| Verschiebungsrichtung | X·a(X) schiebt nach rechts | VL5 Z. 89–91 |
| Welches g bei zwei Möglichkeiten | 1 + X + X³ oder 1 + X² + X³ | Ü6 6_2: beide, äquivalente Codes |
| −1 in F₃ | als −1 oder 2 schreiben | Ü6 6_6: G mit −1; VL5 G₁₁ mit −1 |

#### 14. `kombinatorik.py` – Zählformeln

- **Eingaben:** Modus und Zahlen (n, k; m × n-Gitter; n für Derangements; n für Partitionen).
- **Ausgabe:** Formel und Rechnung: C(n, k), geordnete Auswahl, Kombinationen mit Wiederholung, Gitterwege C(m+n, n), dₙ, p(n).
- **Testfälle:** Ü3 A4 (70), Ü3 A5; Ü4 A2 (d₄ = 9 und Liste); VL2 p(5) = 7 (Z. 340).
- **Konventionen:** keine wesentlichen.

#### 15. `teilbarkeit.py` – Teilbarkeitsregeln und Teilerfunktionen

- **Eingaben:** Zahl x; Modus (9, 11, 7 mit 1001; d(n), σ(n), vollkommen).
- **Ausgabe:** Quersumme, alternierende Quersumme, Dreierblöcke; Formeln für d(n) und σ(n).
- **Testfälle:** Ü1 A1, Ü1 A2 (299792; 10¹⁹ + 1); Ü4 A3–A6 (6, 28, 496; 2¹¹ − 1 = 23·89).
- **Konventionen:** Vorzeichen der alternierenden Quersumme beginnt bei der Einerstelle mit + (Ü1 A1: −2 + 9 − 9 + 7 − 9 + 2).

---

## 5. Fehler und Unsicherheiten in den Quellen

### 5.1 Alle `% UNSICHER`-Kommentare

| Datei, Zeile | Inhalt | Wichtig für die Klausur? |
|---|---|---|
| `uebungen/90_gedaechtnisprotokoll.tex` Z. 46 | GP A3: Werte von a und b fehlen. Nur die Form a·x − b·y = 30 ist bekannt. | Ja. Aufgabentyp üben, Zahlen selbst wählen. Wichtig: rechte Seite 30 ≠ ggT, also hochmultiplizieren |
| `uebungen/90_gedaechtnisprotokoll.tex` Z. 59 | GP A4: Die Wörter „115000711“ und „115000733“ haben nur 9 Stellen. | Ja. Sehr wahrscheinlich fehlt eine 0. Vergleich mit K A5: r₁ = 1150000511, r₂ = 1150000733. Das zweite Wort ist dann identisch |
| `uebungen/90_gedaechtnisprotokoll.tex` Z. 66 | GP A4: Teilaufgabe zweimal „c)“. | Gering. Gemeint ist e) |
| `uebungen/05_uebung.tex` Z. 192 | Ü5 5_5 d: Für r₂ steht x = 13. Richtig ist x = 11 (Spalte (1 2 0)ᵀ). | Mittel. Falscher Testwert für `linearer_code.py`. Das Ergebnis c₂ = 0…0222 stimmt |
| `uebungen/05_uebung.tex` Z. 280 | Ü5 5_8: Im letzten Satz steht [10,8,3]₁₁ statt [10,6,5]₁₁. | Gering. Beide Codes erfüllen Singleton mit Gleichheit |
| `uebungen/06_uebung.tex` Z. 62 | Ü6 6_2: Letzte Zeile von G im 2. Fall ist 0001101. Richtig ist 0001011. | Gering (zyklische Codes). Test muss die korrigierte Zeile nutzen |
| `uebungen/03_uebung.tex` Z. 78 | Ü3 A2: „φ(1000)“ statt „φ(1001)“. | Gering. Nur Schreibfehler, Wert 720 stimmt |
| `wissen/03_potenzen_zm.tex` Z. 254 | VL3: „Alice sendet c an Bob und rechnet“ – gemeint: Bob entschlüsselt. | Mittel. K A1 c fragt: „Wie entschlüsselt **Bob**?“ Antwort: w = cᵈ mod m |
| `wissen/04_lineare_codes.tex` Z. 63 | VL4: „y“ statt r = 010001. | Gering |
| `wissen/04_lineare_codes.tex` Z. 85 | VL4: „w = y“ statt v = w (Abstandsregel). | Gering |
| `wissen/05_zyklische_codes.tex` Z. 334 | VL5: „X⁴ + X³ + X⁴ + X + 1“ statt X⁴ + X³ + X² + X + 1. | Gering |
| `wissen/06_decodierung.tex` Z. 199 | VL6: doppelte Klammer „s(X))“. | Kein Einfluss |

### 5.2 Weitere Auffälligkeiten ohne `UNSICHER`-Markierung (klausurrelevant)

| Stelle | Problem | Folge |
|---|---|---|
| K A1 | m = 43·57, aber 57 = 3·19. m ist kein Produkt zweier Primzahlen. | φ(2451) = 42·2·18 = 1512, d = 821. Die Formel (p−1)(q−1) gibt falsch 2352 und d = 149. Deshalb fordert GP A1 eine „Primzahlzerlegung“ |
| GP A1 | „d = 113“, aber Ü4 A9 mit denselben Zahlen hat e = 113. | Wahrscheinlich falsch erinnert. Wenn doch d = 113: e = 1463 und c = 265 |
| K A5 d, VL6 Z. 279 | Die Vorlesung nennt nur den Fall „s₁ = 0, s₂ ≠ 0“ als nicht decodierbar. K r₁ hat S = (3, 0), also x = 0/3 = 0. | x = 0 ist keine Stelle (Stellen 1 … 10). Also mindestens 2 Fehler. Das muss man selbst begründen |
| K A5 | Die Vorgabe von H ist nicht eindeutig. Grundform und systematische Form geben andere Syndrome. | In der Klausur die Grundform aus VL6/Ü6 6_9 nehmen und dazuschreiben, welche H benutzt wird |
| KV A7 Z. 198, 201 | Schreibweise „d = e^(φ(φ(m))−1) ≡ 1 (mod φ(m))“ ist ungenau. | Gemeint: d ≡ e^(φ(φ(m))−1) mod φ(m) und e·d ≡ 1 mod φ(m) |
| KV A5 | Regel „x + y = 2m“ passt nicht zur Tabelle. | Richtige Regel: x + y ≡ r mod (2m − 1) (Ü1 A3, VL1 Z. 159) |
| Ü4 A1 Z. 13 | Aᵢ ist falsch definiert („i ist Teiler von 1000“). | Richtig: Aᵢ = {k ∈ Ω : i teilt k}. In der Klausur korrekt aufschreiben |
| Ü4 A1 Z. 36 | Produktformel „nur näherungsweise“. | Exakt, wenn N durch alle Teiler-Produkte teilbar ist (K A2: 720, GP A2: 840) |
| VL1 Z. 620 | CRT-Beispiel: „b = 6“ statt b = 1. | Gering; Rechnung nutzt b = 1 |
| Ü5 5_7 | Die Lösung endet ohne Ergebnis. | Hamming-Schranke: erst n = 12 erfüllt 1 + n + n(n−1)/2 ≤ 2^(n−5) (79 ≤ 128) |
| Ü2 A3 Z. 90 | Zweite Gruppentafel mit „ℤ₁₂*“ beschriftet, ist aber ℤ₁₈*. | Gering |
| Ü3 A1 b | Ergebnis nur als −73 mod 60. | Kleinste positive Lösung: 47 |
| H₂(3) in VL4 und Ü5 5_4 | Zwei verschiedene Kontrollmatrizen für denselben Code. | Syndrome hängen von H ab. Immer das H der Aufgabe benutzen |
