# pull the Node.js Docker image
FROM node:alpine

# create the directory inside the container
WORKDIR /app

# copy the package.json files from local machine to the workdir in container
COPY package*.json ./
COPY tsconfig.json ./

# copy the generated modules and all other files to the container
COPY . .

# run npm install in our local machine
RUN npm install
RUN npm run build

# our app is running on port 8080 within the container, so need to expose it
EXPOSE 8080

# the command that starts our app
CMD ["node", "dist/index.js"]