import {defineConfig} from '@playwright/test';
export default defineConfig({testDir:'./tests',testMatch:'*.spec.ts',workers:1,timeout:60000,use:{launchOptions:{executablePath:process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE},baseURL:process.env.WORKBENCH_URL??'http://127.0.0.1:8765',viewport:{width:1500,height:1020},screenshot:'only-on-failure',trace:'retain-on-failure'},reporter:'list'});
