"use strict";

const express = require("express");
const path = require("path");

const app = express();
const port = 3000;
const projectFolder = __dirname;

// Phục vụ index.html, app.js và style.css cho trình duyệt.
app.use(express.static(projectFolder));

app.get("/", function (_request, response) {
  response.sendFile(path.join(projectFolder, "index.html"));
});

app.listen(port, "127.0.0.1", function () {
  console.log("Trang LDR đang chạy tại http://localhost:" + port);
  console.log("Mở đường dẫn trên bằng Chrome hoặc Edge.");
});
