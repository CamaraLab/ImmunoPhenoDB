FROM python:3.8-slim
RUN mkdir /app
WORKDIR /app

# Spark dependencies
# Default values can be overridden at build time
# (ARGS are in lowercase to distinguish them from ENV)

# Use latest spark version (3.5.1 as of 2/29/2024)
ARG spark_version="3.5.1"
ARG hadoop_version="3"
ARG scala_version
# ARG spark_checksum="8883c67e0a138069e597f3e7d4edbbd5c3a565d50b28644aad02856a1ec1da7cb92b8f80454ca427118f69459ea326eaa073cf7b1a860c3b796f4b07c2101319"
ARG openjdk_version="17"

# Install curl and ps
RUN apt-get update --yes && \ 
    apt-get install --yes curl
 
RUN apt-get update --yes && \ 
    apt-get install --yes procps

ENV APACHE_SPARK_VERSION="${spark_version}" \
    HADOOP_VERSION="${hadoop_version}"

# Java Installation
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    "openjdk-${openjdk_version}-jre-headless" \
    ca-certificates-java && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Spark installation
# You need to use https://archive.apache.org/dist/ website if you want to download old Spark versions
# But it seems to be slower, that's why we use the recommended site for download
RUN curl --progress-bar --location --output "spark.tgz" \
        "https://dlcdn.apache.org/spark/spark-${APACHE_SPARK_VERSION}/spark-${APACHE_SPARK_VERSION}-bin-hadoop${HADOOP_VERSION}.tgz" && \
  tar -xf "spark.tgz" -C /usr/local --owner root --group root --no-same-owner && \
  rm "spark.tgz"

# Configure Spark
ENV SPARK_HOME=/usr/local/spark
ENV SPARK_OPTS="--driver-java-options=-Xms1024M --driver-java-options=-Xmx4096M --driver-java-options=-Dlog4j.logLevel=info" \
    PATH="${PATH}:${SPARK_HOME}/bin"

RUN if [ -z "${scala_version}" ]; then \
    ln -s "spark-${APACHE_SPARK_VERSION}-bin-hadoop${HADOOP_VERSION}" "${SPARK_HOME}"; \
  else \
    ln -s "spark-${APACHE_SPARK_VERSION}-bin-hadoop${HADOOP_VERSION}-scala${scala_version}" "${SPARK_HOME}"; \
  fi 

# Not needed, we're not making a jupyter notebook
# Add a link in the before_notebook hook in order to source automatically PYTHONPATH && \
#   ln -s "${SPARK_HOME}/sbin/spark-config.sh" /usr/local/bin/before-notebook.d/10spark-config.sh

ADD requirements.txt /app
RUN pip3 install -r requirements.txt
ADD . /app
EXPOSE 5000
# RUN chmod +x ./entrypoint.sh 
# ENTRYPOINT ["sh", "entrypoint.sh"]
CMD ["gunicorn", "--conf", "/app/gunicorn_config.py", "--bind", "0.0.0.0:5000", "app.wsgi:app"]
