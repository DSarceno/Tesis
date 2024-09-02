# Luego de crear este Dockerfile, desde el directorio donde se encuentra
# ejecutar el siguiente comando:

# docker build -t <nombre de la imagen> .
# donde -t es equivalente a --tag

# Luego se ejecuta el siguiente comando dependiendo del sistema operativo
# Linux/macOS: docker run --rm -v $(pwd):/data compilador-latex
# Windows(PowerShell): docker run --rm -v ${PWD}:/data compilador-latex
# Windows(CMD): docker run --rm -v %cd%:/data compilador-latex




# Se puede abrir docker desde la terminal
# Windows (PowerShell): start-process "C:\Program Files\Docker\Docker\Docker Desktop.exe"
# Sin abrir la interfaz gráfica: docker ps 
# este ultimo mostrara los contenedores en ejecucion


# Para cerrar un contenedor
# docker stop <nombre/ID del contenedor>
# docker stop $(docker ps -q)



# Utiliza una imagen base de LaTeX
FROM blang/latex:ctanfull

# Actualiza el repositorio e instala paquetes adicionales si es necesario
RUN apt-get update && \
    apt-get install -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


WORKDIR /data

CMD pdflatex main.tex && \
    bibtex main.aux && \
    pdflatex main.tex && \
    pdflatex main.tex