# Old dockerfile (for local deployment)
# FROM nginx
# COPY default.conf /etc/nginx/conf.d/default.conf

### NOTE: Update 01-23-2026, the default.conf file bundled with the NGINX docker image is no longer used
### Rather, an updated "default.conf" located on the EC2 instance home directory is used instead to enable HTTPS
### This prevents the need to re-create and re-push a new image each time a modification is made to the conf file

# New dockerfile 04-17-2025 -- bundling React with NGINX on EC2. Must be located in the ROOT directory of project to access the frontend files
# NOTE: Use this dockerfile for creating a NGINX docker image, which will now contain the react front end => to be pushed onto Dockerhub

# Step 1: Use Node compatible version with React
FROM node:18-alpine as build-stage
WORKDIR /app

# Copy package.json and lock file for better caching
# COPY ../vite-project/package.json ../vite-project/package-lock.json* ./
# Try building from root-app
COPY vite-project/package.json vite-project/package-lock.json* ./ 

RUN npm install

# Copy the rest of the React source code
# COPY ../vite-project/ ./
# Try building from root-app
COPY vite-project/ ./

# Run the build script defined in package.json
RUN npm run build
# This should create a build or dist folder inside /app

# Step 2: Serve with Nginx
FROM nginx:1.25-alpine

# Define where NGINX services files from
WORKDIR /usr/share/nginx/html

# Copy static assets from the build stage
# Make sure 'build' matches the output folder name from npm run build
# For vite, I think it's dist (not build)
COPY --from=build-stage /app/dist .

# Copy the custom NGINX configuration (edited here)
# COPY default.conf /etc/nginx/conf.d/default.conf
# Try building from root-app
COPY nginx/default.conf /etc/nginx/conf.d/default.conf

# Expose port 80
EXPOSE 80

# Start NGINX in foreground
CMD ["nginx", "-g", "daemon off;"]