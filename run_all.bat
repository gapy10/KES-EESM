@echo off
REM ============================================================
REM  KES v7 - zagon vseh treh korakov (brez AI)
REM  Faza 1+2: analitika + GA + Pareto + izbor 5 resitev
REM  Faza 3+4: FEMM simulacije (prosti tek + navor)
REM  Zbirni:   Pareto + primerjalni grafi
REM ============================================================

setlocal

REM Premik v mapo, kjer lezi ta skripta (= koren projekta)
cd /d "%~dp0"

echo.
echo === KES v7: zagon vseh treh korakov ===
echo.

REM Izberi Python: poskusi "py", sicer "python"
where py >nul 2>nul
if %ERRORLEVEL%==0 (
    set "PY=py -3"
) else (
    set "PY=python"
)

echo [1/3] Faze 1+2 - analitika + GA + Pareto + izbor 5 resitev...
%PY% -m src.main --config inputs.yaml
if %ERRORLEVEL% neq 0 goto :error

echo.
echo [2/3] Faze 3+4 - FEMM simulacije (prosti tek + navor), lahko traja ~18 min...
echo       Opomba: morebitne vrstice "[info] ... ni zmrezljiv" so NORMALNE -
echo       skript samodejno uporabi najblizjo zmrezljivo resitev s Pareto fronte.
%PY% -m src.run_femm_pareto --step 5
if %ERRORLEVEL% neq 0 goto :error

echo.
echo [3/3] Zbirni Pareto + primerjalni grafi...
%PY% -m src.plot_pareto_femm
if %ERRORLEVEL% neq 0 goto :error

echo.
echo === KONCANO. Rezultati so v mapi outputs\ ===
goto :end

:error
echo.
echo !!! NAPAKA pri zadnjem koraku (izhodna koda %ERRORLEVEL%). Glej sporocila zgoraj. !!!

:end
echo.
pause
endlocal
