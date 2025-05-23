### Passos

0.  Assegurar-se de tenir Python instal·lat fent en una terminal `python --version`.
1. Instal·lar miniconda ([link](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions))
   <details>
   <summary>Procés per Windows</summary>
   1. Obre el CMD i enganxa-hi el següent:
   
      ```
      curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o .\miniconda.exe
      start /wait "" .\miniconda.exe /S
      del .\miniconda.exe
      ```
      Això instal·la l'executable, l'executa i després el suprimeix.
   </details>
2. Afegir miniconda al PATH
   <details>
   <summary>Procés Windows</summary>
   Obre Powe:rShell i fes:
      
      ```
      [System.Environment]::SetEnvironmentVariable("Path", $env:Path + ";$env:USERPROFILE\Miniconda3;$env:USERPROFILE\Miniconda3\Scripts", [System.EnvironmentVariableTarget]::User)
      ```
   Comprova que tot va bé fent:
      ```
      conda --version
      ```

   </details>
   <details>
   <summary>Procés MacOS</summary>

   1. Obre la terminal i fes:
      ```
      nano ~/.zshrc
      ```
   2. Afageix al final de tot una línia que sigui:
      ```
      export PATH="$HOME/Miniconda3/bin:$PATH"
      ```
      Guarda i tanca el fitxer.
   3. Recarrega la configuració fent (a la terminal):
      ```
      source ~/.zshrc
      ```
   4. Comprova que tot va bé fent
      ```
      conda --version
      ```

   </details>
   <details>
   <summary>Procés Linux</summary>

   5. Obre la terminal i fes:
      ```
      nano ~/.bashrc
      ```
   6. Afageix al final de tot una línia que sigui:
      ```
      export PATH="$HOME/Miniconda3/bin:$PATH"
      ```
      Guarda i tanca el fitxer.
   7. Recarrega la configuració fent (a la terminal):
      ```
      source ~/.bashrc
      ```
   8. Comprova que tot va bé fent
      ```
      conda --version
      ```

   </details>
3. Reiniciar la terminal
4. Fer `conda init` i després `conda activate base`
   
   Amb això estem obrint l'entorn base (també podríem crear un entorn personalitzat per cada projecte però ho farem en el base per simplicitat). Allà hi instal·larem tots els paquets i llibreries.
5. Instal·lar els següents paquets:
   1. pandas
   2. scipy
   3. numpy
   4. matplotlib
   5. cartopy
   6. xarray
   7. netcdf4
   
   Els instal·lem utilitzant: `conda install [nomPaquet]`.
6. Baixar aquest repositori com a .zip, descompirmir-lo, guardar-lo a algun lloc (Documents per exemple) i obrir-lo en el  VSCode.
7. Obrir el fitxer `scripts/spatial_calculation_onemember_numpy.ipynb`.
8. A dalt a la dreta del Jupyter, clicar on posa "Select Kernel" i triar "Python Environments" > "Base (Python 3.13.2)". 
      
   Nota: El número de la versió dependrà de quina versió de Python tingueu instal·lada. Això indica que el Kernel de Jupyter s'executarà en l'entorn 'base' de miniconda (que és on hem instal·lat tots els paquets).
9.  Executar la 1a cel·la (tardarà una estona perquè haurà d'instal·larà el kernel de Python i algunes extensions com Jupyter).
10. Si la 1a cel·la s'executa correctament la part difícil ja està feta. Això vol dir que estan totes les llibreries instal·lades correctament.
11. Anar al Drive ([link](https://drive.google.com/drive/folders/1XAAp8DirojrlLGHW-w2-fX8Ob4S9AK4C)) i baixar l'arxiu corresponent al membre que us ha tocat. Posar-lo dins la carpeta `member_data`. 
    
    Nota: L'arxiu que ja hi ha dins la carpeta és a mode d'exemple, és un fitxer buit. Aquests fitxers de dades pesen massa per poder-los pujar a GitHub, però en local no hi ha problemes.
12. A la cel·la 2 canviar la línia
    ```py
    numMembre='18' 
    ```
    Pel número que us ha tocat.
13. Executa la cel·la 2 i seguidament la 3. Si no dona error ja estàs llegint correctament el fitxer de dades.
14. Executa les següents cel·les, en la 7 s'hauria de generar un mapa de la climatologia.
15. Si fins aquí no hi ha hagut cap error, ja tot va bé. Executa la resta de cel·les una a una.
