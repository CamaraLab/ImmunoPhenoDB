import express, { Express, Request, Response } from 'express'
import dotenv from 'dotenv';
import router from './routes/db_routes.js'

dotenv.config({ path: '../.env'})

import cors from 'cors';

const app: Express = express();
app.use(express.json());
app.use(cors());

const NODE_PORT = 8080; // You can also use process.env.NODE_PORT if needed
const SERVER_TIMEOUT = 20 * 60 * 1000; // 20 minutes in milliseconds

// Start the server and capture the server instance
const server = app.listen(NODE_PORT, () => {
    console.log(`Node server is up and running at port ${NODE_PORT}`);
});

// Set the server timeout to 20 minutes
server.timeout = SERVER_TIMEOUT;

app.use('/api', router);