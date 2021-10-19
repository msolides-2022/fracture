docker rm -f notebook && \
docker run --name notebook -w /home/fenics -v $(pwd):/home/fenics/shared -d -p 127.0.0.1:8888:8888 cmaurini/fenics-stable-20.04 'jupyter-notebook --ip=0.0.0.0 --allow-root' 
docker logs notebook
