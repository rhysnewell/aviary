import os
import pysam
import numpy
import glob


def write_css(outfile):
    with open(outfile, 'w') as o:
        o.write('''@import url(https://fonts.googleapis.com/css?family=Poppins:300,400,700);
:root {
  --blue: #007bff;
  --indigo: #6610f2;
  --purple: #6f42c1;
  --pink: #e83e8c;
  --red: #dc3545;
  --orange: #fd7e14;
  --yellow: #ffc107;
  --green: #28a745;
  --teal: #20c997;
  --cyan: #17a2b8;
  --white: #fff;
  --gray: #6c757d;
  --gray-dark: #343a40;
  --primary: #512479;
  --secondary: #333;
  --success: #28a745;
  --info: #17a2b8;
  --warning: #ffc107;
  --danger: #dc3545;
  --light: #f8f9fa;
  --dark: #343a40;
  --breakpoint-xs: 0;
  --breakpoint-sm: 576px;
  --breakpoint-md: 768px;
  --breakpoint-lg: 992px;
  --breakpoint-xl: 1200px;
  --font-family-sans-serif: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";
  --font-family-monospace: SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; }
*,
*::before,
*::after {
  box-sizing: border-box; }
.toggle-vis {
  color: green !important;
  cursor: pointer;
}
.greened {
  color: red !important;
}
html {
  font-family: sans-serif;
  line-height: 1.15;
  -webkit-text-size-adjust: 100%;
  -ms-text-size-adjust: 100%;
  -ms-overflow-style: scrollbar;
  -webkit-tap-highlight-color: rgba(0, 0, 0, 0); }
@-ms-viewport {
  width: device-width; }
article, aside, dialog, figcaption, figure, footer, header, hgroup, main, nav, section {
  display: block; }
body {
  margin: 0;
  font-family: "Poppins", sans-serif;
  font-size: 1rem;
  font-weight: 300;
  line-height: 2.4;
  color: #999;
  text-align: left;
  background-color: #fff; }
[tabindex="-1"]:focus {
  outline: 0 !important; }
hr {
  box-sizing: content-box;
  height: 0;
  overflow: visible; }
h1, h2, h3, h4, h5, h6 {
  margin-top: 0;
  margin-bottom: 0.5rem; }
p {
  margin-top: 0;
  margin-bottom: 1rem; }
abbr[title],
abbr[data-original-title] {
  text-decoration: underline;
  text-decoration: underline dotted;
  cursor: help;
  border-bottom: 0; }
address {
  margin-bottom: 1rem;
  font-style: normal;
  line-height: inherit; }
ol,
ul,
dl {
  margin-top: 0;
  margin-bottom: 1rem; }
ol ol,
ul ul,
ol ul,
ul ol {
  margin-bottom: 0; }
dt {
  font-weight: 700; }
dd {
  margin-bottom: .5rem;
  margin-left: 0; }
blockquote {
  margin: 0 0 1rem; }
dfn {
  font-style: italic; }
b,
strong {
  font-weight: bolder; }
small {
  font-size: 80%; }
sub,
sup {
  position: relative;
  font-size: 75%;
  line-height: 0;
  vertical-align: baseline; }
sub {
  bottom: -.25em; }
sup {
  top: -.5em; }
a {
  color: #512479;
  text-decoration: none;
  background-color: transparent;
  -webkit-text-decoration-skip: objects; }
  a:hover {
    color: #523047;
    text-decoration: underline; }
a:not([href]):not([tabindex]) {
  color: inherit;
  text-decoration: none; }
  a:not([href]):not([tabindex]):hover, a:not([href]):not([tabindex]):focus {
    color: inherit;
    text-decoration: none; }
  a:not([href]):not([tabindex]):focus {
    outline: 0; }
pre,
code,
kbd,
samp {
  font-family: monospace, monospace;
  font-size: 1em; }
pre {
  margin-top: 0;
  margin-bottom: 1rem;
  overflow: auto;
  -ms-overflow-style: scrollbar; }
figure {
  margin: 0 0 1rem; }
img {
  vertical-align: middle;
  border-style: none; }
svg:not(:root) {
  overflow: hidden; }
table {
  border-collapse: collapse; }
caption {
  padding-top: 0.75rem;
  padding-bottom: 0.75rem;
  color: #6c757d;
  text-align: left;
  caption-side: bottom; }
th {
  text-align: inherit; }
label {
  display: inline-block;
  margin-bottom: .5rem; }
button {
  border-radius: 0; }
button:focus {
  outline: 1px dotted;
  outline: 5px auto -webkit-focus-ring-color; }
input,
button,
select,
optgroup,
textarea {
  margin: 0;
  font-family: inherit;
  font-size: inherit;
  line-height: inherit; }
button,
input {
  overflow: visible; }
button,
select {
  text-transform: none; }
button,
html [type="button"],
[type="reset"],
[type="submit"] {
  -webkit-appearance: button; }
button::-moz-focus-inner,
[type="button"]::-moz-focus-inner,
[type="reset"]::-moz-focus-inner,
[type="submit"]::-moz-focus-inner {
  padding: 0;
  border-style: none; }
input[type="radio"],
input[type="checkbox"] {
  box-sizing: border-box;
  padding: 0; }
input[type="date"],
input[type="time"],
input[type="datetime-local"],
input[type="month"] {
  -webkit-appearance: listbox; }
textarea {
  overflow: auto;
  resize: vertical; }
fieldset {
  min-width: 0;
  padding: 0;
  margin: 0;
  border: 0; }
legend {
  display: block;
  width: 100%;
  max-width: 100%;
  padding: 0;
  margin-bottom: .5rem;
  font-size: 1.5rem;
  line-height: inherit;
  color: inherit;
  white-space: normal; }
progress {
  vertical-align: baseline; }
[type="number"]::-webkit-inner-spin-button,
[type="number"]::-webkit-outer-spin-button {
  height: auto; }
[type="search"] {
  outline-offset: -2px;
  -webkit-appearance: none; }
[type="search"]::-webkit-search-cancel-button,
[type="search"]::-webkit-search-decoration {
  -webkit-appearance: none; }
::-webkit-file-upload-button {
  font: inherit;
  -webkit-appearance: button; }
output {
  display: inline-block; }
summary {
  display: list-item;
  cursor: pointer; }
template {
  display: none; }
[hidden] {
  display: none !important; }
h1, h2, h3, h4, h5, h6,
.h1, .h2, .h3, .h4, .h5, .h6 {
  margin-bottom: 0.5rem;
  font-family: inherit;
  font-weight: 500;
  line-height: 1.2;
  color: inherit; }
h1, .h1 {
  font-size: 2.5rem; }
h2, .h2 {
  font-size: 2rem; }
h3, .h3 {
  font-size: 1.75rem; }
h4, .h4 {
  font-size: 1.5rem; }
h5, .h5 {
  font-size: 1.25rem; }
h6, .h6 {
  font-size: 1rem; }
.lead {
  font-size: 1.25rem;
  font-weight: 300; }
.display-1 {
  font-size: 6rem;
  font-weight: 300;
  line-height: 1.2; }
.display-2 {
  font-size: 5.5rem;
  font-weight: 300;
  line-height: 1.2; }
.display-3 {
  font-size: 4.5rem;
  font-weight: 300;
  line-height: 1.2; }
.display-4 {
  font-size: 3.5rem;
  font-weight: 300;
  line-height: 1.2; }
hr {
  margin-top: 1rem;
  margin-bottom: 1rem;
  border: 0;
  border-top: 1px solid rgba(0, 0, 0, 0.1); }
small,
.small {
  font-size: 80%;
  font-weight: 400; }
mark,
.mark {
  padding: 0.2em;
  background-color: #fcf8e3; }
.list-unstyled {
  padding-left: 0;
  list-style: none; }
.list-inline {
  padding-left: 0;
  list-style: none; }
.list-inline-item {
  display: inline-block; }
  .list-inline-item:not(:last-child) {
    margin-right: 0.5rem; }
.initialism {
  font-size: 90%;
  text-transform: uppercase; }
.blockquote {
  margin-bottom: 1rem;
  font-size: 1.25rem; }
.blockquote-footer {
  display: block;
  font-size: 80%;
  color: #6c757d; }
  .blockquote-footer::before {
    content: "\2014 \00A0"; }
.img-fluid {
  max-width: 100%;
  height: auto; }
.img-thumbnail {
  padding: 0.25rem;
  background-color: #fff;
  border: 1px solid #dee2e6;
  border-radius: 0px;
  max-width: 100%;
  height: auto; }
.figure {
  display: inline-block; }
.figure-img {
  margin-bottom: 0.5rem;
  line-height: 1; }
.figure-caption {
  font-size: 90%;
  color: #6c757d; }
code,
kbd,
pre,
samp {
  font-family: SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; }
code {
  font-size: 87.5%;
  color: #e83e8c;
  word-break: break-word; }
  a > code {
    color: inherit; }
kbd {
  padding: 0.2rem 0.4rem;
  font-size: 87.5%;
  color: #fff;
  background-color: #212529;
  border-radius: 0px; }
  kbd kbd {
    padding: 0;
    font-size: 100%;
    font-weight: 700; }
pre {
  display: block;
  font-size: 87.5%;
  color: #212529; }
  pre code {
    font-size: inherit;
    color: inherit;
    word-break: normal; }
.pre-scrollable {
  max-height: 340px;
  overflow-y: scroll; }
.container {
  width: 100%;
  padding-right: 15px;
  padding-left: 15px;
  margin-right: auto;
  margin-left: auto; }
  @media (min-width: 576px) {
    .container {
      max-width: 740px; } }
  @media (min-width: 768px) {
    .container {
      max-width: 960px; } }
  @media (min-width: 992px) {
    .container {
      max-width: 1060px; } }
  @media (min-width: 1200px) {
    .container {
      max-width: 1280px; } }
.container-fluid {
  width: 100%;
  padding-right: 15px;
  padding-left: 15px;
  margin-right: auto;
  margin-left: auto; }
.row {
  display: flex;
  flex-wrap: wrap;
  margin-right: -15px;
  margin-left: -15px; }
.no-gutters {
  margin-right: 0;
  margin-left: 0; }
  .no-gutters > .col,
  .no-gutters > [class*="col-"] {
    padding-right: 0;
    padding-left: 0; }
.col-1, .col-2, .col-3, .col-4, .col-5, .col-6, .col-7, .col-8, .col-9, .col-10, .col-11, .col-12, .col,
.col-auto, .col-sm-1, .col-sm-2, .col-sm-3, .col-sm-4, .col-sm-5, .col-sm-6, .col-sm-7, .col-sm-8, .col-sm-9, .col-sm-10, .col-sm-11, .col-sm-12, .col-sm,
.col-sm-auto, .col-md-1, .col-md-2, .col-md-3, .col-md-4, .col-md-5, .col-md-6, .col-md-7, .col-md-8, .col-md-9, .col-md-10, .col-md-11, .col-md-12, .col-md,
.col-md-auto, .col-lg-1, .col-lg-2, .col-lg-3, .col-lg-4, .col-lg-5, .col-lg-6, .col-lg-7, .col-lg-8, .col-lg-9, .col-lg-10, .col-lg-11, .col-lg-12, .col-lg,
.col-lg-auto, .col-xl-1, .col-xl-2, .col-xl-3, .col-xl-4, .col-xl-5, .col-xl-6, .col-xl-7, .col-xl-8, .col-xl-9, .col-xl-10, .col-xl-11, .col-xl-12, .col-xl,
.col-xl-auto {
  position: relative;
  width: 100%;
  min-height: 1px;
  padding-right: 15px;
  padding-left: 15px; }
.col {
  flex-basis: 0;
  flex-grow: 1;
  max-width: 100%; }
.col-auto {
  flex: 0 0 auto;
  width: auto;
  max-width: none; }
.col-1 {
  flex: 0 0 8.3333333333%;
  max-width: 8.3333333333%; }
.col-2 {
  flex: 0 0 16.6666666667%;
  max-width: 16.6666666667%; }
.col-3 {
  flex: 0 0 25%;
  max-width: 25%; }
.col-4 {
  flex: 0 0 33.3333333333%;
  max-width: 33.3333333333%; }
.col-5 {
  flex: 0 0 41.6666666667%;
  max-width: 41.6666666667%; }
.col-6 {
  flex: 0 0 50%;
  max-width: 50%; }
.col-7 {
  flex: 0 0 58.3333333333%;
  max-width: 58.3333333333%; }
.col-8 {
  flex: 0 0 66.6666666667%;
  max-width: 66.6666666667%; }
.col-9 {
  flex: 0 0 75%;
  max-width: 75%; }
.col-10 {
  flex: 0 0 83.3333333333%;
  max-width: 83.3333333333%; }
.col-11 {
  flex: 0 0 91.6666666667%;
  max-width: 91.6666666667%; }
.col-12 {
  flex: 0 0 100%;
  max-width: 100%; }
.order-first {
  order: -1; }
.order-last {
  order: 13; }
.order-0 {
  order: 0; }
.order-1 {
  order: 1; }
.order-2 {
  order: 2; }
.order-3 {
  order: 3; }
.order-4 {
  order: 4; }
.order-5 {
  order: 5; }
.order-6 {
  order: 6; }
.order-7 {
  order: 7; }
.order-8 {
  order: 8; }
.order-9 {
  order: 9; }
.order-10 {
  order: 10; }
.order-11 {
  order: 11; }
.order-12 {
  order: 12; }
.offset-1 {
  margin-left: 8.3333333333%; }
.offset-2 {
  margin-left: 16.6666666667%; }
.offset-3 {
  margin-left: 25%; }
.offset-4 {
  margin-left: 33.3333333333%; }
.offset-5 {
  margin-left: 41.6666666667%; }
.offset-6 {
  margin-left: 50%; }
.offset-7 {
  margin-left: 58.3333333333%; }
.offset-8 {
  margin-left: 66.6666666667%; }
.offset-9 {
  margin-left: 75%; }
.offset-10 {
  margin-left: 83.3333333333%; }
.offset-11 {
  margin-left: 91.6666666667%; }
@media (min-width: 576px) {
  .col-sm {
    flex-basis: 0;
    flex-grow: 1;
    max-width: 100%; }
  .col-sm-auto {
    flex: 0 0 auto;
    width: auto;
    max-width: none; }
  .col-sm-1 {
    flex: 0 0 8.3333333333%;
    max-width: 8.3333333333%; }
  .col-sm-2 {
    flex: 0 0 16.6666666667%;
    max-width: 16.6666666667%; }
  .col-sm-3 {
    flex: 0 0 25%;
    max-width: 25%; }
  .col-sm-4 {
    flex: 0 0 33.3333333333%;
    max-width: 33.3333333333%; }
  .col-sm-5 {
    flex: 0 0 41.6666666667%;
    max-width: 41.6666666667%; }
  .col-sm-6 {
    flex: 0 0 50%;
    max-width: 50%; }
  .col-sm-7 {
    flex: 0 0 58.3333333333%;
    max-width: 58.3333333333%; }
  .col-sm-8 {
    flex: 0 0 66.6666666667%;
    max-width: 66.6666666667%; }
  .col-sm-9 {
    flex: 0 0 75%;
    max-width: 75%; }
  .col-sm-10 {
    flex: 0 0 83.3333333333%;
    max-width: 83.3333333333%; }
  .col-sm-11 {
    flex: 0 0 91.6666666667%;
    max-width: 91.6666666667%; }
  .col-sm-12 {
    flex: 0 0 100%;
    max-width: 100%; }
  .order-sm-first {
    order: -1; }
  .order-sm-last {
    order: 13; }
  .order-sm-0 {
    order: 0; }
  .order-sm-1 {
    order: 1; }
  .order-sm-2 {
    order: 2; }
  .order-sm-3 {
    order: 3; }
  .order-sm-4 {
    order: 4; }
  .order-sm-5 {
    order: 5; }
  .order-sm-6 {
    order: 6; }
  .order-sm-7 {
    order: 7; }
  .order-sm-8 {
    order: 8; }
  .order-sm-9 {
    order: 9; }
  .order-sm-10 {
    order: 10; }
  .order-sm-11 {
    order: 11; }
  .order-sm-12 {
    order: 12; }
  .offset-sm-0 {
    margin-left: 0; }
  .offset-sm-1 {
    margin-left: 8.3333333333%; }
  .offset-sm-2 {
    margin-left: 16.6666666667%; }
  .offset-sm-3 {
    margin-left: 25%; }
  .offset-sm-4 {
    margin-left: 33.3333333333%; }
  .offset-sm-5 {
    margin-left: 41.6666666667%; }
  .offset-sm-6 {
    margin-left: 50%; }
  .offset-sm-7 {
    margin-left: 58.3333333333%; }
  .offset-sm-8 {
    margin-left: 66.6666666667%; }
  .offset-sm-9 {
    margin-left: 75%; }
  .offset-sm-10 {
    margin-left: 83.3333333333%; }
  .offset-sm-11 {
    margin-left: 91.6666666667%; } }
@media (min-width: 768px) {
  .col-md {
    flex-basis: 0;
    flex-grow: 1;
    max-width: 100%; }
  .col-md-auto {
    flex: 0 0 auto;
    width: auto;
    max-width: none; }
  .col-md-1 {
    flex: 0 0 8.3333333333%;
    max-width: 8.3333333333%; }
  .col-md-2 {
    flex: 0 0 16.6666666667%;
    max-width: 16.6666666667%; }
  .col-md-3 {
    flex: 0 0 25%;
    max-width: 25%; }
  .col-md-4 {
    flex: 0 0 33.3333333333%;
    max-width: 33.3333333333%; }
  .col-md-5 {
    flex: 0 0 41.6666666667%;
    max-width: 41.6666666667%; }
  .col-md-6 {
    flex: 0 0 50%;
    max-width: 50%; }
  .col-md-7 {
    flex: 0 0 58.3333333333%;
    max-width: 58.3333333333%; }
  .col-md-8 {
    flex: 0 0 66.6666666667%;
    max-width: 66.6666666667%; }
  .col-md-9 {
    flex: 0 0 75%;
    max-width: 75%; }
  .col-md-10 {
    flex: 0 0 83.3333333333%;
    max-width: 83.3333333333%; }
  .col-md-11 {
    flex: 0 0 91.6666666667%;
    max-width: 91.6666666667%; }
  .col-md-12 {
    flex: 0 0 100%;
    max-width: 100%; }
  .order-md-first {
    order: -1; }
  .order-md-last {
    order: 13; }
  .order-md-0 {
    order: 0; }
  .order-md-1 {
    order: 1; }
  .order-md-2 {
    order: 2; }
  .order-md-3 {
    order: 3; }
  .order-md-4 {
    order: 4; }
  .order-md-5 {
    order: 5; }
  .order-md-6 {
    order: 6; }
  .order-md-7 {
    order: 7; }
  .order-md-8 {
    order: 8; }
  .order-md-9 {
    order: 9; }
  .order-md-10 {
    order: 10; }
  .order-md-11 {
    order: 11; }
  .order-md-12 {
    order: 12; }
  .offset-md-0 {
    margin-left: 0; }
  .offset-md-1 {
    margin-left: 8.3333333333%; }
  .offset-md-2 {
    margin-left: 16.6666666667%; }
  .offset-md-3 {
    margin-left: 25%; }
  .offset-md-4 {
    margin-left: 33.3333333333%; }
  .offset-md-5 {
    margin-left: 41.6666666667%; }
  .offset-md-6 {
    margin-left: 50%; }
  .offset-md-7 {
    margin-left: 58.3333333333%; }
  .offset-md-8 {
    margin-left: 66.6666666667%; }
  .offset-md-9 {
    margin-left: 75%; }
  .offset-md-10 {
    margin-left: 83.3333333333%; }
  .offset-md-11 {
    margin-left: 91.6666666667%; } }
@media (min-width: 992px) {
  .col-lg {
    flex-basis: 0;
    flex-grow: 1;
    max-width: 100%; }
  .col-lg-auto {
    flex: 0 0 auto;
    width: auto;
    max-width: none; }
  .col-lg-1 {
    flex: 0 0 8.3333333333%;
    max-width: 8.3333333333%; }
  .col-lg-2 {
    flex: 0 0 16.6666666667%;
    max-width: 16.6666666667%; }
  .col-lg-3 {
    flex: 0 0 25%;
    max-width: 25%; }
  .col-lg-4 {
    flex: 0 0 33.3333333333%;
    max-width: 33.3333333333%; }
  .col-lg-5 {
    flex: 0 0 41.6666666667%;
    max-width: 41.6666666667%; }
  .col-lg-6 {
    flex: 0 0 50%;
    max-width: 50%; }
  .col-lg-7 {
    flex: 0 0 58.3333333333%;
    max-width: 58.3333333333%; }
  .col-lg-8 {
    flex: 0 0 66.6666666667%;
    max-width: 66.6666666667%; }
  .col-lg-9 {
    flex: 0 0 75%;
    max-width: 75%; }
  .col-lg-10 {
    flex: 0 0 83.3333333333%;
    max-width: 83.3333333333%; }
  .col-lg-11 {
    flex: 0 0 91.6666666667%;
    max-width: 91.6666666667%; }
  .col-lg-12 {
    flex: 0 0 100%;
    max-width: 100%; }
  .order-lg-first {
    order: -1; }
  .order-lg-last {
    order: 13; }
  .order-lg-0 {
    order: 0; }
  .order-lg-1 {
    order: 1; }
  .order-lg-2 {
    order: 2; }
  .order-lg-3 {
    order: 3; }
  .order-lg-4 {
    order: 4; }
  .order-lg-5 {
    order: 5; }
  .order-lg-6 {
    order: 6; }
  .order-lg-7 {
    order: 7; }
  .order-lg-8 {
    order: 8; }
  .order-lg-9 {
    order: 9; }
  .order-lg-10 {
    order: 10; }
  .order-lg-11 {
    order: 11; }
  .order-lg-12 {
    order: 12; }
  .offset-lg-0 {
    margin-left: 0; }
  .offset-lg-1 {
    margin-left: 8.3333333333%; }
  .offset-lg-2 {
    margin-left: 16.6666666667%; }
  .offset-lg-3 {
    margin-left: 25%; }
  .offset-lg-4 {
    margin-left: 33.3333333333%; }
  .offset-lg-5 {
    margin-left: 41.6666666667%; }
  .offset-lg-6 {
    margin-left: 50%; }
  .offset-lg-7 {
    margin-left: 58.3333333333%; }
  .offset-lg-8 {
    margin-left: 66.6666666667%; }
  .offset-lg-9 {
    margin-left: 75%; }
  .offset-lg-10 {
    margin-left: 83.3333333333%; }
  .offset-lg-11 {
    margin-left: 91.6666666667%; } }
@media (min-width: 1200px) {
  .col-xl {
    flex-basis: 0;
    flex-grow: 1;
    max-width: 100%; }
  .col-xl-auto {
    flex: 0 0 auto;
    width: auto;
    max-width: none; }
  .col-xl-1 {
    flex: 0 0 8.3333333333%;
    max-width: 8.3333333333%; }
  .col-xl-2 {
    flex: 0 0 16.6666666667%;
    max-width: 16.6666666667%; }
  .col-xl-3 {
    flex: 0 0 25%;
    max-width: 25%; }
  .col-xl-4 {
    flex: 0 0 33.3333333333%;
    max-width: 33.3333333333%; }
  .col-xl-5 {
    flex: 0 0 41.6666666667%;
    max-width: 41.6666666667%; }
  .col-xl-6 {
    flex: 0 0 50%;
    max-width: 50%; }
  .col-xl-7 {
    flex: 0 0 58.3333333333%;
    max-width: 58.3333333333%; }
  .col-xl-8 {
    flex: 0 0 66.6666666667%;
    max-width: 66.6666666667%; }
  .col-xl-9 {
    flex: 0 0 75%;
    max-width: 75%; }
  .col-xl-10 {
    flex: 0 0 83.3333333333%;
    max-width: 83.3333333333%; }
  .col-xl-11 {
    flex: 0 0 91.6666666667%;
    max-width: 91.6666666667%; }
  .col-xl-12 {
    flex: 0 0 100%;
    max-width: 100%; }
  .order-xl-first {
    order: -1; }
  .order-xl-last {
    order: 13; }
  .order-xl-0 {
    order: 0; }
  .order-xl-1 {
    order: 1; }
  .order-xl-2 {
    order: 2; }
  .order-xl-3 {
    order: 3; }
  .order-xl-4 {
    order: 4; }
  .order-xl-5 {
    order: 5; }
  .order-xl-6 {
    order: 6; }
  .order-xl-7 {
    order: 7; }
  .order-xl-8 {
    order: 8; }
  .order-xl-9 {
    order: 9; }
  .order-xl-10 {
    order: 10; }
  .order-xl-11 {
    order: 11; }
  .order-xl-12 {
    order: 12; }
  .offset-xl-0 {
    margin-left: 0; }
  .offset-xl-1 {
    margin-left: 8.3333333333%; }
  .offset-xl-2 {
    margin-left: 16.6666666667%; }
  .offset-xl-3 {
    margin-left: 25%; }
  .offset-xl-4 {
    margin-left: 33.3333333333%; }
  .offset-xl-5 {
    margin-left: 41.6666666667%; }
  .offset-xl-6 {
    margin-left: 50%; }
  .offset-xl-7 {
    margin-left: 58.3333333333%; }
  .offset-xl-8 {
    margin-left: 66.6666666667%; }
  .offset-xl-9 {
    margin-left: 75%; }
  .offset-xl-10 {
    margin-left: 83.3333333333%; }
  .offset-xl-11 {
    margin-left: 91.6666666667%; } }
.table {
  width: 100%;
  max-width: 100%;
  margin-bottom: 1rem;
  font-size: 14px;
  color: #000000;
  background-color: transparent; }
  .table th,
  .table td {
    padding: 0.75rem;
    vertical-align: top;
    border-top: 1px solid #dee2e6; }
  .table thead th {
    vertical-align: bottom;
    border-bottom: 2px solid #dee2e6; }
  .table tbody + tbody {
    border-top: 2px solid #dee2e6; }
  .table .table {
    background-color: #fff; }
tr:hover {background-color: #f0f0f0;}
@media (max-width: 575.98px) {
  .table-responsive-sm {
    display: block;
    width: 100%;
    overflow-x: auto;
    -webkit-overflow-scrolling: touch;
    -ms-overflow-style: -ms-autohiding-scrollbar; }
    .table-responsive-sm > .table-bordered {
      border: 0; } }
@media (max-width: 767.98px) {
  .table-responsive-md {
    display: block;
    width: 100%;
    overflow-x: auto;
    -webkit-overflow-scrolling: touch;
    -ms-overflow-style: -ms-autohiding-scrollbar; }
    .table-responsive-md > .table-bordered {
      border: 0; } }
@media (max-width: 991.98px) {
  .table-responsive-lg {
    display: block;
    width: 100%;
    overflow-x: auto;
    -webkit-overflow-scrolling: touch;
    -ms-overflow-style: -ms-autohiding-scrollbar; }
    .table-responsive-lg > .table-bordered {
      border: 0; } }
@media (max-width: 1199.98px) {
  .table-responsive-xl {
    display: block;
    width: 100%;
    overflow-x: auto;
    -webkit-overflow-scrolling: touch;
    -ms-overflow-style: -ms-autohiding-scrollbar; }
    .table-responsive-xl > .table-bordered {
      border: 0; } }
.table-responsive {
  display: block;
  width: 100%;
  overflow-x: auto;
  -webkit-overflow-scrolling: touch;
  -ms-overflow-style: -ms-autohiding-scrollbar; }
  .table-responsive > .table-bordered {
    border: 0; }
.form-control {
  display: block;
  width: 100%;
  padding: 0.375rem 0.75rem;
  font-size: 1rem;
  line-height: 2.4;
  color: #495057;
  background-color: #fff;
  background-clip: padding-box;
  border: 1px solid #ced4da;
  border-radius: 0px;
  transition: border-color 0.15s ease-in-out, box-shadow 0.15s ease-in-out; }
  .form-control::-ms-expand {
    background-color: transparent;
    border: 0; }
  .form-control:focus {
    color: #495057;
    background-color: #fff;
    border-color: #737373;
    outline: 0;
    box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.25); }
  .form-control::placeholder {
    color: #6c757d;
    opacity: 1; }
  .form-control:disabled, .form-control[readonly] {
    background-color: #e9ecef;
    opacity: 1; }
select.form-control:not([size]):not([multiple]) {
  height: calc(3.15rem + 2px); }
select.form-control:focus::-ms-value {
  color: #495057;
  background-color: #fff; }
.form-control-file,
.form-control-range {
  display: block;
  width: 100%; }
.col-form-label {
  padding-top: calc(0.375rem + 1px);
  padding-bottom: calc(0.375rem + 1px);
  margin-bottom: 0;
  font-size: inherit;
  line-height: 2.4; }
.col-form-label-lg {
  padding-top: calc(0.5rem + 1px);
  padding-bottom: calc(0.5rem + 1px);
  font-size: 1.25rem;
  line-height: 1.5; }
.col-form-label-sm {
  padding-top: calc(0.25rem + 1px);
  padding-bottom: calc(0.25rem + 1px);
  font-size: 0.875rem;
  line-height: 1.5; }
.form-control-plaintext {
  display: block;
  width: 100%;
  padding-top: 0.375rem;
  padding-bottom: 0.375rem;
  margin-bottom: 0;
  line-height: 2.4;
  background-color: transparent;
  border: solid transparent;
  border-width: 1px 0; }
  .form-control-plaintext.form-control-sm, .input-group-sm > .form-control-plaintext.form-control,
  .input-group-sm > .input-group-prepend > .form-control-plaintext.input-group-text,
  .input-group-sm > .input-group-append > .form-control-plaintext.input-group-text,
  .input-group-sm > .input-group-prepend > .form-control-plaintext.btn,
  .input-group-sm > .input-group-append > .form-control-plaintext.btn, .form-control-plaintext.form-control-lg, .input-group-lg > .form-control-plaintext.form-control,
  .input-group-lg > .input-group-prepend > .form-control-plaintext.input-group-text,
  .input-group-lg > .input-group-append > .form-control-plaintext.input-group-text,
  .input-group-lg > .input-group-prepend > .form-control-plaintext.btn,
  .input-group-lg > .input-group-append > .form-control-plaintext.btn {
    padding-right: 0;
    padding-left: 0; }
.form-control-sm, .input-group-sm > .form-control,
.input-group-sm > .input-group-prepend > .input-group-text,
.input-group-sm > .input-group-append > .input-group-text,
.input-group-sm > .input-group-prepend > .btn,
.input-group-sm > .input-group-append > .btn {
  padding: 0.25rem 0.5rem;
  font-size: 0.875rem;
  line-height: 1.5;
  border-radius: 0px; }
select.form-control-sm:not([size]):not([multiple]), .input-group-sm > select.form-control:not([size]):not([multiple]),
.input-group-sm > .input-group-prepend > select.input-group-text:not([size]):not([multiple]),
.input-group-sm > .input-group-append > select.input-group-text:not([size]):not([multiple]),
.input-group-sm > .input-group-prepend > select.btn:not([size]):not([multiple]),
.input-group-sm > .input-group-append > select.btn:not([size]):not([multiple]) {
  height: calc(1.8125rem + 2px); }
.form-control-lg, .input-group-lg > .form-control,
.input-group-lg > .input-group-prepend > .input-group-text,
.input-group-lg > .input-group-append > .input-group-text,
.input-group-lg > .input-group-prepend > .btn,
.input-group-lg > .input-group-append > .btn {
  padding: 0.5rem 1rem;
  font-size: 1.25rem;
  line-height: 1.5;
  border-radius: 0px; }
select.form-control-lg:not([size]):not([multiple]), .input-group-lg > select.form-control:not([size]):not([multiple]),
.input-group-lg > .input-group-prepend > select.input-group-text:not([size]):not([multiple]),
.input-group-lg > .input-group-append > select.input-group-text:not([size]):not([multiple]),
.input-group-lg > .input-group-prepend > select.btn:not([size]):not([multiple]),
.input-group-lg > .input-group-append > select.btn:not([size]):not([multiple]) {
  height: calc(2.875rem + 2px); }
.form-group {
  margin-bottom: 1rem; }
.form-text {
  display: block;
  margin-top: 0.25rem; }
.form-row {
  display: flex;
  flex-wrap: wrap;
  margin-right: -5px;
  margin-left: -5px; }
  .form-row > .col,
  .form-row > [class*="col-"] {
    padding-right: 5px;
    padding-left: 5px; }
.form-check {
  position: relative;
  display: block;
  padding-left: 1.25rem; }
.form-check-input {
  position: absolute;
  margin-top: 0.3rem;
  margin-left: -1.25rem; }
  .form-check-input:disabled ~ .form-check-label {
    color: #6c757d; }
.form-check-label {
  margin-bottom: 0; }
.form-check-inline {
  display: inline-flex;
  align-items: center;
  padding-left: 0;
  margin-right: 0.75rem; }
  .form-check-inline .form-check-input {
    position: static;
    margin-top: 0;
    margin-right: 0.3125rem;
    margin-left: 0; }
.valid-feedback {
  display: none;
  width: 100%;
  margin-top: 0.25rem;
  font-size: 80%;
  color: #28a745; }
.valid-tooltip {
  position: absolute;
  top: 100%;
  z-index: 5;
  display: none;
  max-width: 100%;
  padding: .5rem;
  margin-top: .1rem;
  font-size: .875rem;
  line-height: 1;
  color: #fff;
  background-color: rgba(40, 167, 69, 0.8);
  border-radius: .2rem; }
.was-validated .form-control:valid, .form-control.is-valid,
.was-validated .custom-select:valid,
.custom-select.is-valid {
  border-color: #28a745; }
  .was-validated .form-control:valid:focus, .form-control.is-valid:focus,
  .was-validated .custom-select:valid:focus,
  .custom-select.is-valid:focus {
    border-color: #28a745;
    box-shadow: 0 0 0 0.2rem rgba(40, 167, 69, 0.25); }
  .was-validated .form-control:valid ~ .valid-feedback,
  .was-validated .form-control:valid ~ .valid-tooltip, .form-control.is-valid ~ .valid-feedback,
  .form-control.is-valid ~ .valid-tooltip,
  .was-validated .custom-select:valid ~ .valid-feedback,
  .was-validated .custom-select:valid ~ .valid-tooltip,
  .custom-select.is-valid ~ .valid-feedback,
  .custom-select.is-valid ~ .valid-tooltip {
    display: block; }
.was-validated .form-check-input:valid ~ .form-check-label, .form-check-input.is-valid ~ .form-check-label {
  color: #28a745; }
.was-validated .form-check-input:valid ~ .valid-feedback,
.was-validated .form-check-input:valid ~ .valid-tooltip, .form-check-input.is-valid ~ .valid-feedback,
.form-check-input.is-valid ~ .valid-tooltip {
  display: block; }
.was-validated .custom-control-input:valid ~ .custom-control-label, .custom-control-input.is-valid ~ .custom-control-label {
  color: #28a745; }
  .was-validated .custom-control-input:valid ~ .custom-control-label::before, .custom-control-input.is-valid ~ .custom-control-label::before {
    background-color: #71dd8a; }
.was-validated .custom-control-input:valid ~ .valid-feedback,
.was-validated .custom-control-input:valid ~ .valid-tooltip, .custom-control-input.is-valid ~ .valid-feedback,
.custom-control-input.is-valid ~ .valid-tooltip {
  display: block; }
.was-validated .custom-control-input:valid:checked ~ .custom-control-label::before, .custom-control-input.is-valid:checked ~ .custom-control-label::before {
  background-color: #34ce57; }
.was-validated .custom-control-input:valid:focus ~ .custom-control-label::before, .custom-control-input.is-valid:focus ~ .custom-control-label::before {
  box-shadow: 0 0 0 1px #fff, 0 0 0 0.2rem rgba(40, 167, 69, 0.25); }
.was-validated .custom-file-input:valid ~ .custom-file-label, .custom-file-input.is-valid ~ .custom-file-label {
  border-color: #28a745; }
  .was-validated .custom-file-input:valid ~ .custom-file-label::before, .custom-file-input.is-valid ~ .custom-file-label::before {
    border-color: inherit; }
.was-validated .custom-file-input:valid ~ .valid-feedback,
.was-validated .custom-file-input:valid ~ .valid-tooltip, .custom-file-input.is-valid ~ .valid-feedback,
.custom-file-input.is-valid ~ .valid-tooltip {
  display: block; }
.was-validated .custom-file-input:valid:focus ~ .custom-file-label, .custom-file-input.is-valid:focus ~ .custom-file-label {
  box-shadow: 0 0 0 0.2rem rgba(40, 167, 69, 0.25); }
.invalid-feedback {
  display: none;
  width: 100%;
  margin-top: 0.25rem;
  font-size: 80%;
  color: #dc3545; }
.invalid-tooltip {
  position: absolute;
  top: 100%;
  z-index: 5;
  display: none;
  max-width: 100%;
  padding: .5rem;
  margin-top: .1rem;
  font-size: .875rem;
  line-height: 1;
  color: #fff;
  background-color: rgba(220, 53, 69, 0.8);
  border-radius: .2rem; }
.was-validated .form-control:invalid, .form-control.is-invalid,
.was-validated .custom-select:invalid,
.custom-select.is-invalid {
  border-color: #dc3545; }
  .was-validated .form-control:invalid:focus, .form-control.is-invalid:focus,
  .was-validated .custom-select:invalid:focus,
  .custom-select.is-invalid:focus {
    border-color: #dc3545;
    box-shadow: 0 0 0 0.2rem rgba(220, 53, 69, 0.25); }
  .was-validated .form-control:invalid ~ .invalid-feedback,
  .was-validated .form-control:invalid ~ .invalid-tooltip, .form-control.is-invalid ~ .invalid-feedback,
  .form-control.is-invalid ~ .invalid-tooltip,
  .was-validated .custom-select:invalid ~ .invalid-feedback,
  .was-validated .custom-select:invalid ~ .invalid-tooltip,
  .custom-select.is-invalid ~ .invalid-feedback,
  .custom-select.is-invalid ~ .invalid-tooltip {
    display: block; }
.was-validated .form-check-input:invalid ~ .form-check-label, .form-check-input.is-invalid ~ .form-check-label {
  color: #dc3545; }
.was-validated .form-check-input:invalid ~ .invalid-feedback,
.was-validated .form-check-input:invalid ~ .invalid-tooltip, .form-check-input.is-invalid ~ .invalid-feedback,
.form-check-input.is-invalid ~ .invalid-tooltip {
  display: block; }
.was-validated .custom-control-input:invalid ~ .custom-control-label, .custom-control-input.is-invalid ~ .custom-control-label {
  color: #dc3545; }
  .was-validated .custom-control-input:invalid ~ .custom-control-label::before, .custom-control-input.is-invalid ~ .custom-control-label::before {
    background-color: #efa2a9; }
.was-validated .custom-control-input:invalid ~ .invalid-feedback,
.was-validated .custom-control-input:invalid ~ .invalid-tooltip, .custom-control-input.is-invalid ~ .invalid-feedback,
.custom-control-input.is-invalid ~ .invalid-tooltip {
  display: block; }
.was-validated .custom-control-input:invalid:checked ~ .custom-control-label::before, .custom-control-input.is-invalid:checked ~ .custom-control-label::before {
  background-color: #e4606d; }
.was-validated .custom-control-input:invalid:focus ~ .custom-control-label::before, .custom-control-input.is-invalid:focus ~ .custom-control-label::before {
  box-shadow: 0 0 0 1px #fff, 0 0 0 0.2rem rgba(220, 53, 69, 0.25); }
.was-validated .custom-file-input:invalid ~ .custom-file-label, .custom-file-input.is-invalid ~ .custom-file-label {
  border-color: #dc3545; }
  .was-validated .custom-file-input:invalid ~ .custom-file-label::before, .custom-file-input.is-invalid ~ .custom-file-label::before {
    border-color: inherit; }
.was-validated .custom-file-input:invalid ~ .invalid-feedback,
.was-validated .custom-file-input:invalid ~ .invalid-tooltip, .custom-file-input.is-invalid ~ .invalid-feedback,
.custom-file-input.is-invalid ~ .invalid-tooltip {
  display: block; }
.was-validated .custom-file-input:invalid:focus ~ .custom-file-label, .custom-file-input.is-invalid:focus ~ .custom-file-label {
  box-shadow: 0 0 0 0.2rem rgba(220, 53, 69, 0.25); }
.form-inline {
  display: flex;
  flex-flow: row wrap;
  align-items: center; }
  .form-inline .form-check {
    width: 100%; }
  @media (min-width: 576px) {
    .form-inline label {
      display: flex;
      align-items: center;
      justify-content: center;
      margin-bottom: 0; }
    .form-inline .form-group {
      display: flex;
      flex: 0 0 auto;
      flex-flow: row wrap;
      align-items: center;
      margin-bottom: 0; }
    .form-inline .form-control {
      display: inline-block;
      width: auto;
      vertical-align: middle; }
    .form-inline .form-control-plaintext {
      display: inline-block; }
    .form-inline .input-group {
      width: auto; }
    .form-inline .form-check {
      display: flex;
      align-items: center;
      justify-content: center;
      width: auto;
      padding-left: 0; }
    .form-inline .form-check-input {
      position: relative;
      margin-top: 0;
      margin-right: 0.25rem;
      margin-left: 0; }
    .form-inline .custom-control {
      align-items: center;
      justify-content: center; }
    .form-inline .custom-control-label {
      margin-bottom: 0; } }
.btn {
  display: inline-block;
  font-weight: 400;
  text-align: center;
  white-space: nowrap;
  vertical-align: middle;
  user-select: none;
  border: 1px solid transparent;
  padding: 0.375rem 0.75rem;
  font-size: 1rem;
  line-height: 2.4;
  border-radius: 0px;
  transition: color 0.15s ease-in-out, background-color 0.15s ease-in-out, border-color 0.15s ease-in-out, box-shadow 0.15s ease-in-out; }
  .btn:hover, .btn:focus {
    text-decoration: none; }
  .btn:focus, .btn.focus {
    outline: 0;
    box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.25); }
  .btn.disabled, .btn:disabled {
    opacity: 0.65; }
  .btn:not(:disabled):not(.disabled) {
    cursor: pointer; }
  .btn:not(:disabled):not(.disabled):active, .btn:not(:disabled):not(.disabled).active {
    background-image: none; }
a.btn.disabled,
fieldset:disabled a.btn {
  pointer-events: none; }
.btn-primary {
  color: #fff;
  background-color: #512479;
  border-color: #512479; }
  .btn-primary:hover {
    color: #fff;
    background-color: #6a3e5c;
    border-color: #49075e; }
  .btn-primary:focus, .btn-primary.focus {
    box-shadow: 0 0 0 0.2rem rgba(130, 76, 113, 0.5); }
  .btn-primary.disabled, .btn-primary:disabled {
    color: #fff;
    background-color: #512479;
    border-color: #512479; }
  .btn-primary:not(:disabled):not(.disabled):active, .btn-primary:not(:disabled):not(.disabled).active, .show > .btn-primary.dropdown-toggle {
    color: #fff;
    background-color: #49075e;
    border-color: #5a344e; }
    .btn-primary:not(:disabled):not(.disabled):active:focus, .btn-primary:not(:disabled):not(.disabled).active:focus, .show > .btn-primary.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(130, 76, 113, 0.5); }
.btn-secondary {
  color: #fff;
  background-color: #333;
  border-color: #333; }
  .btn-secondary:hover {
    color: #fff;
    background-color: #202020;
    border-color: #1a1a1a; }
  .btn-secondary:focus, .btn-secondary.focus {
    box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.5); }
  .btn-secondary.disabled, .btn-secondary:disabled {
    color: #fff;
    background-color: #333;
    border-color: #333; }
  .btn-secondary:not(:disabled):not(.disabled):active, .btn-secondary:not(:disabled):not(.disabled).active, .show > .btn-secondary.dropdown-toggle {
    color: #fff;
    background-color: #1a1a1a;
    border-color: #131313; }
    .btn-secondary:not(:disabled):not(.disabled):active:focus, .btn-secondary:not(:disabled):not(.disabled).active:focus, .show > .btn-secondary.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.5); }
.btn-success {
  color: #fff;
  background-color: #28a745;
  border-color: #28a745; }
  .btn-success:hover {
    color: #fff;
    background-color: #218838;
    border-color: #1e7e34; }
  .btn-success:focus, .btn-success.focus {
    box-shadow: 0 0 0 0.2rem rgba(40, 167, 69, 0.5); }
  .btn-success.disabled, .btn-success:disabled {
    color: #fff;
    background-color: #28a745;
    border-color: #28a745; }
  .btn-success:not(:disabled):not(.disabled):active, .btn-success:not(:disabled):not(.disabled).active, .show > .btn-success.dropdown-toggle {
    color: #fff;
    background-color: #1e7e34;
    border-color: #1c7430; }
    .btn-success:not(:disabled):not(.disabled):active:focus, .btn-success:not(:disabled):not(.disabled).active:focus, .show > .btn-success.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(40, 167, 69, 0.5); }
.btn-info {
  color: #fff;
  background-color: #17a2b8;
  border-color: #17a2b8; }
  .btn-info:hover {
    color: #fff;
    background-color: #138496;
    border-color: #117a8b; }
  .btn-info:focus, .btn-info.focus {
    box-shadow: 0 0 0 0.2rem rgba(23, 162, 184, 0.5); }
  .btn-info.disabled, .btn-info:disabled {
    color: #fff;
    background-color: #17a2b8;
    border-color: #17a2b8; }
  .btn-info:not(:disabled):not(.disabled):active, .btn-info:not(:disabled):not(.disabled).active, .show > .btn-info.dropdown-toggle {
    color: #fff;
    background-color: #117a8b;
    border-color: #10707f; }
    .btn-info:not(:disabled):not(.disabled):active:focus, .btn-info:not(:disabled):not(.disabled).active:focus, .show > .btn-info.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(23, 162, 184, 0.5); }
.btn-warning {
  color: #212529;
  background-color: #ffc107;
  border-color: #ffc107; }
  .btn-warning:hover {
    color: #212529;
    background-color: #e0a800;
    border-color: #d39e00; }
  .btn-warning:focus, .btn-warning.focus {
    box-shadow: 0 0 0 0.2rem rgba(255, 193, 7, 0.5); }
  .btn-warning.disabled, .btn-warning:disabled {
    color: #212529;
    background-color: #ffc107;
    border-color: #ffc107; }
  .btn-warning:not(:disabled):not(.disabled):active, .btn-warning:not(:disabled):not(.disabled).active, .show > .btn-warning.dropdown-toggle {
    color: #212529;
    background-color: #d39e00;
    border-color: #c69500; }
    .btn-warning:not(:disabled):not(.disabled):active:focus, .btn-warning:not(:disabled):not(.disabled).active:focus, .show > .btn-warning.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(255, 193, 7, 0.5); }
.btn-danger {
  color: #fff;
  background-color: #dc3545;
  border-color: #dc3545; }
  .btn-danger:hover {
    color: #fff;
    background-color: #c82333;
    border-color: #bd2130; }
  .btn-danger:focus, .btn-danger.focus {
    box-shadow: 0 0 0 0.2rem rgba(220, 53, 69, 0.5); }
  .btn-danger.disabled, .btn-danger:disabled {
    color: #fff;
    background-color: #dc3545;
    border-color: #dc3545; }
  .btn-danger:not(:disabled):not(.disabled):active, .btn-danger:not(:disabled):not(.disabled).active, .show > .btn-danger.dropdown-toggle {
    color: #fff;
    background-color: #bd2130;
    border-color: #b21f2d; }
    .btn-danger:not(:disabled):not(.disabled):active:focus, .btn-danger:not(:disabled):not(.disabled).active:focus, .show > .btn-danger.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(220, 53, 69, 0.5); }
.btn-light {
  color: #212529;
  background-color: #f8f9fa;
  border-color: #f8f9fa; }
  .btn-light:hover {
    color: #212529;
    background-color: #e2e6ea;
    border-color: #dae0e5; }
  .btn-light:focus, .btn-light.focus {
    box-shadow: 0 0 0 0.2rem rgba(248, 249, 250, 0.5); }
  .btn-light.disabled, .btn-light:disabled {
    color: #212529;
    background-color: #f8f9fa;
    border-color: #f8f9fa; }
  .btn-light:not(:disabled):not(.disabled):active, .btn-light:not(:disabled):not(.disabled).active, .show > .btn-light.dropdown-toggle {
    color: #212529;
    background-color: #dae0e5;
    border-color: #d3d9df; }
    .btn-light:not(:disabled):not(.disabled):active:focus, .btn-light:not(:disabled):not(.disabled).active:focus, .show > .btn-light.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(248, 249, 250, 0.5); }
.btn-dark {
  color: #fff;
  background-color: #343a40;
  border-color: #343a40; }
  .btn-dark:hover {
    color: #fff;
    background-color: #23272b;
    border-color: #1d2124; }
  .btn-dark:focus, .btn-dark.focus {
    box-shadow: 0 0 0 0.2rem rgba(52, 58, 64, 0.5); }
  .btn-dark.disabled, .btn-dark:disabled {
    color: #fff;
    background-color: #343a40;
    border-color: #343a40; }
  .btn-dark:not(:disabled):not(.disabled):active, .btn-dark:not(:disabled):not(.disabled).active, .show > .btn-dark.dropdown-toggle {
    color: #fff;
    background-color: #1d2124;
    border-color: #171a1d; }
    .btn-dark:not(:disabled):not(.disabled):active:focus, .btn-dark:not(:disabled):not(.disabled).active:focus, .show > .btn-dark.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(52, 58, 64, 0.5); }
.btn-outline-primary {
  color: #512479;
  background-color: transparent;
  background-image: none;
  border-color: #512479; }
  .btn-outline-primary:hover {
    color: #fff;
    background-color: #512479;
    border-color: #512479; }
  .btn-outline-primary:focus, .btn-outline-primary.focus {
    box-shadow: 0 0 0 0.2rem rgba(130, 76, 113, 0.5); }
  .btn-outline-primary.disabled, .btn-outline-primary:disabled {
    color: #512479;
    background-color: transparent; }
  .btn-outline-primary:not(:disabled):not(.disabled):active, .btn-outline-primary:not(:disabled):not(.disabled).active, .show > .btn-outline-primary.dropdown-toggle {
    color: #fff;
    background-color: #512479;
    border-color: #512479; }
    .btn-outline-primary:not(:disabled):not(.disabled):active:focus, .btn-outline-primary:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-primary.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(130, 76, 113, 0.5); }
.btn-outline-secondary {
  color: #333;
  background-color: transparent;
  background-image: none;
  border-color: #333; }
  .btn-outline-secondary:hover {
    color: #fff;
    background-color: #333;
    border-color: #333; }
  .btn-outline-secondary:focus, .btn-outline-secondary.focus {
    box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.5); }
  .btn-outline-secondary.disabled, .btn-outline-secondary:disabled {
    color: #333;
    background-color: transparent; }
  .btn-outline-secondary:not(:disabled):not(.disabled):active, .btn-outline-secondary:not(:disabled):not(.disabled).active, .show > .btn-outline-secondary.dropdown-toggle {
    color: #fff;
    background-color: #333;
    border-color: #333; }
    .btn-outline-secondary:not(:disabled):not(.disabled):active:focus, .btn-outline-secondary:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-secondary.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.5); }
.btn-outline-success {
  color: #28a745;
  background-color: transparent;
  background-image: none;
  border-color: #28a745; }
  .btn-outline-success:hover {
    color: #fff;
    background-color: #28a745;
    border-color: #28a745; }
  .btn-outline-success:focus, .btn-outline-success.focus {
    box-shadow: 0 0 0 0.2rem rgba(40, 167, 69, 0.5); }
  .btn-outline-success.disabled, .btn-outline-success:disabled {
    color: #28a745;
    background-color: transparent; }
  .btn-outline-success:not(:disabled):not(.disabled):active, .btn-outline-success:not(:disabled):not(.disabled).active, .show > .btn-outline-success.dropdown-toggle {
    color: #fff;
    background-color: #28a745;
    border-color: #28a745; }
    .btn-outline-success:not(:disabled):not(.disabled):active:focus, .btn-outline-success:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-success.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(40, 167, 69, 0.5); }
.btn-outline-info {
  color: #17a2b8;
  background-color: transparent;
  background-image: none;
  border-color: #17a2b8; }
  .btn-outline-info:hover {
    color: #fff;
    background-color: #17a2b8;
    border-color: #17a2b8; }
  .btn-outline-info:focus, .btn-outline-info.focus {
    box-shadow: 0 0 0 0.2rem rgba(23, 162, 184, 0.5); }
  .btn-outline-info.disabled, .btn-outline-info:disabled {
    color: #17a2b8;
    background-color: transparent; }
  .btn-outline-info:not(:disabled):not(.disabled):active, .btn-outline-info:not(:disabled):not(.disabled).active, .show > .btn-outline-info.dropdown-toggle {
    color: #fff;
    background-color: #17a2b8;
    border-color: #17a2b8; }
    .btn-outline-info:not(:disabled):not(.disabled):active:focus, .btn-outline-info:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-info.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(23, 162, 184, 0.5); }
.btn-outline-warning {
  color: #ffc107;
  background-color: transparent;
  background-image: none;
  border-color: #ffc107; }
  .btn-outline-warning:hover {
    color: #212529;
    background-color: #ffc107;
    border-color: #ffc107; }
  .btn-outline-warning:focus, .btn-outline-warning.focus {
    box-shadow: 0 0 0 0.2rem rgba(255, 193, 7, 0.5); }
  .btn-outline-warning.disabled, .btn-outline-warning:disabled {
    color: #ffc107;
    background-color: transparent; }
  .btn-outline-warning:not(:disabled):not(.disabled):active, .btn-outline-warning:not(:disabled):not(.disabled).active, .show > .btn-outline-warning.dropdown-toggle {
    color: #212529;
    background-color: #ffc107;
    border-color: #ffc107; }
    .btn-outline-warning:not(:disabled):not(.disabled):active:focus, .btn-outline-warning:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-warning.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(255, 193, 7, 0.5); }
.btn-outline-danger {
  color: #dc3545;
  background-color: transparent;
  background-image: none;
  border-color: #dc3545; }
  .btn-outline-danger:hover {
    color: #fff;
    background-color: #dc3545;
    border-color: #dc3545; }
  .btn-outline-danger:focus, .btn-outline-danger.focus {
    box-shadow: 0 0 0 0.2rem rgba(220, 53, 69, 0.5); }
  .btn-outline-danger.disabled, .btn-outline-danger:disabled {
    color: #dc3545;
    background-color: transparent; }
  .btn-outline-danger:not(:disabled):not(.disabled):active, .btn-outline-danger:not(:disabled):not(.disabled).active, .show > .btn-outline-danger.dropdown-toggle {
    color: #fff;
    background-color: #dc3545;
    border-color: #dc3545; }
    .btn-outline-danger:not(:disabled):not(.disabled):active:focus, .btn-outline-danger:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-danger.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(220, 53, 69, 0.5); }
.btn-outline-light {
  color: #f8f9fa;
  background-color: transparent;
  background-image: none;
  border-color: #f8f9fa; }
  .btn-outline-light:hover {
    color: #212529;
    background-color: #f8f9fa;
    border-color: #f8f9fa; }
  .btn-outline-light:focus, .btn-outline-light.focus {
    box-shadow: 0 0 0 0.2rem rgba(248, 249, 250, 0.5); }
  .btn-outline-light.disabled, .btn-outline-light:disabled {
    color: #f8f9fa;
    background-color: transparent; }
  .btn-outline-light:not(:disabled):not(.disabled):active, .btn-outline-light:not(:disabled):not(.disabled).active, .show > .btn-outline-light.dropdown-toggle {
    color: #212529;
    background-color: #f8f9fa;
    border-color: #f8f9fa; }
    .btn-outline-light:not(:disabled):not(.disabled):active:focus, .btn-outline-light:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-light.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(248, 249, 250, 0.5); }
.btn-outline-dark {
  color: #343a40;
  background-color: transparent;
  background-image: none;
  border-color: #343a40; }
  .btn-outline-dark:hover {
    color: #fff;
    background-color: #343a40;
    border-color: #343a40; }
  .btn-outline-dark:focus, .btn-outline-dark.focus {
    box-shadow: 0 0 0 0.2rem rgba(52, 58, 64, 0.5); }
  .btn-outline-dark.disabled, .btn-outline-dark:disabled {
    color: #343a40;
    background-color: transparent; }
  .btn-outline-dark:not(:disabled):not(.disabled):active, .btn-outline-dark:not(:disabled):not(.disabled).active, .show > .btn-outline-dark.dropdown-toggle {
    color: #fff;
    background-color: #343a40;
    border-color: #343a40; }
    .btn-outline-dark:not(:disabled):not(.disabled):active:focus, .btn-outline-dark:not(:disabled):not(.disabled).active:focus, .show > .btn-outline-dark.dropdown-toggle:focus {
      box-shadow: 0 0 0 0.2rem rgba(52, 58, 64, 0.5); }
.btn-link {
  font-weight: 400;
  color: #512479;
  background-color: transparent; }
  .btn-link:hover {
    color: #523047;
    text-decoration: underline;
    background-color: transparent;
    border-color: transparent; }
  .btn-link:focus, .btn-link.focus {
    text-decoration: underline;
    border-color: transparent;
    box-shadow: none; }
  .btn-link:disabled, .btn-link.disabled {
    color: #6c757d; }
.btn-lg, .btn-group-lg > .btn {
  padding: 0.5rem 1rem;
  font-size: 1.25rem;
  line-height: 1.5;
  border-radius: 0px; }
.btn-sm, .btn-group-sm > .btn {
  padding: 0.25rem 0.5rem;
  font-size: 0.875rem;
  line-height: 1.5;
  border-radius: 0px; }
.btn-block {
  display: block;
  width: 100%; }
  .btn-block + .btn-block {
    margin-top: 0.5rem; }
input[type="submit"].btn-block,
input[type="reset"].btn-block,
input[type="button"].btn-block {
  width: 100%; }
.fade {
  opacity: 0;
  transition: opacity 0.15s linear; }
  .fade.show {
    opacity: 1; }
.collapse {
  display: none; }
  .collapse.show {
    display: block; }
tr.collapse.show {
  display: table-row; }
tbody.collapse.show {
  display: table-row-group; }
.collapsing {
  position: relative;
  height: 0;
  overflow: hidden;
  transition: height 0.35s ease; }
.dropup,
.dropdown {
  position: relative; }
.dropdown-toggle::after {
  display: inline-block;
  width: 0;
  height: 0;
  margin-left: 0.255em;
  vertical-align: 0.255em;
  content: "";
  border-top: 0.3em solid;
  border-right: 0.3em solid transparent;
  border-bottom: 0;
  border-left: 0.3em solid transparent; }
.dropdown-toggle:empty::after {
  margin-left: 0; }
.dropdown-menu {
  position: absolute;
  top: 100%;
  left: 0;
  z-index: 1000;
  display: none;
  float: left;
  min-width: 10rem;
  padding: 0.5rem 0;
  margin: 0.125rem 0 0;
  font-size: 1rem;
  color: #999;
  text-align: left;
  list-style: none;
  background-color: #fff;
  background-clip: padding-box;
  border: 1px solid rgba(0, 0, 0, 0.15);
  border-radius: 0px; }
.dropup .dropdown-menu {
  margin-top: 0;
  margin-bottom: 0.125rem; }
.dropup .dropdown-toggle::after {
  display: inline-block;
  width: 0;
  height: 0;
  margin-left: 0.255em;
  vertical-align: 0.255em;
  content: "";
  border-top: 0;
  border-right: 0.3em solid transparent;
  border-bottom: 0.3em solid;
  border-left: 0.3em solid transparent; }
.dropup .dropdown-toggle:empty::after {
  margin-left: 0; }
.dropright .dropdown-menu {
  margin-top: 0;
  margin-left: 0.125rem; }
.dropright .dropdown-toggle::after {
  display: inline-block;
  width: 0;
  height: 0;
  margin-left: 0.255em;
  vertical-align: 0.255em;
  content: "";
  border-top: 0.3em solid transparent;
  border-bottom: 0.3em solid transparent;
  border-left: 0.3em solid; }
.dropright .dropdown-toggle:empty::after {
  margin-left: 0; }
.dropright .dropdown-toggle::after {
  vertical-align: 0; }
.dropleft .dropdown-menu {
  margin-top: 0;
  margin-right: 0.125rem; }
.dropleft .dropdown-toggle::after {
  display: inline-block;
  width: 0;
  height: 0;
  margin-left: 0.255em;
  vertical-align: 0.255em;
  content: ""; }
.dropleft .dropdown-toggle::after {
  display: none; }
.dropleft .dropdown-toggle::before {
  display: inline-block;
  width: 0;
  height: 0;
  margin-right: 0.255em;
  vertical-align: 0.255em;
  content: "";
  border-top: 0.3em solid transparent;
  border-right: 0.3em solid;
  border-bottom: 0.3em solid transparent; }
.dropleft .dropdown-toggle:empty::after {
  margin-left: 0; }
.dropleft .dropdown-toggle::before {
  vertical-align: 0; }
.dropdown-divider {
  height: 0;
  margin: 0.5rem 0;
  overflow: hidden;
  border-top: 1px solid #e9ecef; }
.dropdown-item {
  display: block;
  width: 100%;
  padding: 0.25rem 1.5rem;
  clear: both;
  font-weight: 400;
  color: #212529;
  text-align: inherit;
  white-space: nowrap;
  background-color: transparent;
  border: 0; }
  .dropdown-item:hover, .dropdown-item:focus {
    color: #16181b;
    text-decoration: none;
    background-color: #f8f9fa; }
  .dropdown-item.active, .dropdown-item:active {
    color: #fff;
    text-decoration: none;
    background-color: #333; }
  .dropdown-item.disabled, .dropdown-item:disabled {
    color: #6c757d;
    background-color: transparent; }
.dropdown-menu.show {
  display: block; }
.dropdown-header {
  display: block;
  padding: 0.5rem 1.5rem;
  margin-bottom: 0;
  font-size: 0.875rem;
  color: #6c757d;
  white-space: nowrap; }
.btn-group,
.btn-group-vertical {
  position: relative;
  display: inline-flex;
  vertical-align: middle; }
  .btn-group > .btn,
  .btn-group-vertical > .btn {
    position: relative;
    flex: 0 1 auto; }
    .btn-group > .btn:hover,
    .btn-group-vertical > .btn:hover {
      z-index: 1; }
    .btn-group > .btn:focus, .btn-group > .btn:active, .btn-group > .btn.active,
    .btn-group-vertical > .btn:focus,
    .btn-group-vertical > .btn:active,
    .btn-group-vertical > .btn.active {
      z-index: 1; }
  .btn-group .btn + .btn,
  .btn-group .btn + .btn-group,
  .btn-group .btn-group + .btn,
  .btn-group .btn-group + .btn-group,
  .btn-group-vertical .btn + .btn,
  .btn-group-vertical .btn + .btn-group,
  .btn-group-vertical .btn-group + .btn,
  .btn-group-vertical .btn-group + .btn-group {
    margin-left: -1px; }
.btn-toolbar {
  display: flex;
  flex-wrap: wrap;
  justify-content: flex-start; }
  .btn-toolbar .input-group {
    width: auto; }
.btn-group > .btn:first-child {
  margin-left: 0; }
.btn-group > .btn:not(:last-child):not(.dropdown-toggle),
.btn-group > .btn-group:not(:last-child) > .btn {
  border-top-right-radius: 0;
  border-bottom-right-radius: 0; }
.btn-group > .btn:not(:first-child),
.btn-group > .btn-group:not(:first-child) > .btn {
  border-top-left-radius: 0;
  border-bottom-left-radius: 0; }
.dropdown-toggle-split {
  padding-right: 0.5625rem;
  padding-left: 0.5625rem; }
  .dropdown-toggle-split::after {
    margin-left: 0; }
.btn-sm + .dropdown-toggle-split, .btn-group-sm > .btn + .dropdown-toggle-split {
  padding-right: 0.375rem;
  padding-left: 0.375rem; }
.btn-lg + .dropdown-toggle-split, .btn-group-lg > .btn + .dropdown-toggle-split {
  padding-right: 0.75rem;
  padding-left: 0.75rem; }
.btn-group-vertical {
  flex-direction: column;
  align-items: flex-start;
  justify-content: center; }
  .btn-group-vertical .btn,
  .btn-group-vertical .btn-group {
    width: 100%; }
  .btn-group-vertical > .btn + .btn,
  .btn-group-vertical > .btn + .btn-group,
  .btn-group-vertical > .btn-group + .btn,
  .btn-group-vertical > .btn-group + .btn-group {
    margin-top: -1px;
    margin-left: 0; }
  .btn-group-vertical > .btn:not(:last-child):not(.dropdown-toggle),
  .btn-group-vertical > .btn-group:not(:last-child) > .btn {
    border-bottom-right-radius: 0;
    border-bottom-left-radius: 0; }
  .btn-group-vertical > .btn:not(:first-child),
  .btn-group-vertical > .btn-group:not(:first-child) > .btn {
    border-top-left-radius: 0;
    border-top-right-radius: 0; }
.btn-group-toggle > .btn,
.btn-group-toggle > .btn-group > .btn {
  margin-bottom: 0; }
  .btn-group-toggle > .btn input[type="radio"],
  .btn-group-toggle > .btn input[type="checkbox"],
  .btn-group-toggle > .btn-group > .btn input[type="radio"],
  .btn-group-toggle > .btn-group > .btn input[type="checkbox"] {
    position: absolute;
    clip: rect(0, 0, 0, 0);
    pointer-events: none; }
.input-group {
  position: relative;
  display: flex;
  flex-wrap: wrap;
  align-items: stretch;
  width: 100%; }
  .input-group > .form-control,
  .input-group > .custom-select,
  .input-group > .custom-file {
    position: relative;
    flex: 1 1 auto;
    width: 1%;
    margin-bottom: 0; }
    .input-group > .form-control:focus,
    .input-group > .custom-select:focus,
    .input-group > .custom-file:focus {
      z-index: 3; }
    .input-group > .form-control + .form-control,
    .input-group > .form-control + .custom-select,
    .input-group > .form-control + .custom-file,
    .input-group > .custom-select + .form-control,
    .input-group > .custom-select + .custom-select,
    .input-group > .custom-select + .custom-file,
    .input-group > .custom-file + .form-control,
    .input-group > .custom-file + .custom-select,
    .input-group > .custom-file + .custom-file {
      margin-left: -1px; }
  .input-group > .form-control:not(:last-child),
  .input-group > .custom-select:not(:last-child) {
    border-top-right-radius: 0;
    border-bottom-right-radius: 0; }
  .input-group > .form-control:not(:first-child),
  .input-group > .custom-select:not(:first-child) {
    border-top-left-radius: 0;
    border-bottom-left-radius: 0; }
  .input-group > .custom-file {
    display: flex;
    align-items: center; }
    .input-group > .custom-file:not(:last-child) .custom-file-label, .input-group > .custom-file:not(:last-child) .custom-file-label::before {
      border-top-right-radius: 0;
      border-bottom-right-radius: 0; }
    .input-group > .custom-file:not(:first-child) .custom-file-label, .input-group > .custom-file:not(:first-child) .custom-file-label::before {
      border-top-left-radius: 0;
      border-bottom-left-radius: 0; }
.input-group-prepend,
.input-group-append {
  display: flex; }
  .input-group-prepend .btn,
  .input-group-append .btn {
    position: relative;
    z-index: 2; }
  .input-group-prepend .btn + .btn,
  .input-group-prepend .btn + .input-group-text,
  .input-group-prepend .input-group-text + .input-group-text,
  .input-group-prepend .input-group-text + .btn,
  .input-group-append .btn + .btn,
  .input-group-append .btn + .input-group-text,
  .input-group-append .input-group-text + .input-group-text,
  .input-group-append .input-group-text + .btn {
    margin-left: -1px; }
.input-group-prepend {
  margin-right: -1px; }
.input-group-append {
  margin-left: -1px; }
.input-group-text {
  display: flex;
  align-items: center;
  padding: 0.375rem 0.75rem;
  margin-bottom: 0;
  font-size: 1rem;
  font-weight: 400;
  line-height: 2.4;
  color: #495057;
  text-align: center;
  white-space: nowrap;
  background-color: #e9ecef;
  border: 1px solid #ced4da;
  border-radius: 0px; }
  .input-group-text input[type="radio"],
  .input-group-text input[type="checkbox"] {
    margin-top: 0; }
.input-group > .input-group-prepend > .btn,
.input-group > .input-group-prepend > .input-group-text,
.input-group > .input-group-append:not(:last-child) > .btn,
.input-group > .input-group-append:not(:last-child) > .input-group-text,
.input-group > .input-group-append:last-child > .btn:not(:last-child):not(.dropdown-toggle),
.input-group > .input-group-append:last-child > .input-group-text:not(:last-child) {
  border-top-right-radius: 0;
  border-bottom-right-radius: 0; }
.input-group > .input-group-append > .btn,
.input-group > .input-group-append > .input-group-text,
.input-group > .input-group-prepend:not(:first-child) > .btn,
.input-group > .input-group-prepend:not(:first-child) > .input-group-text,
.input-group > .input-group-prepend:first-child > .btn:not(:first-child),
.input-group > .input-group-prepend:first-child > .input-group-text:not(:first-child) {
  border-top-left-radius: 0;
  border-bottom-left-radius: 0; }
.custom-control {
  position: relative;
  display: block;
  min-height: 2.4rem;
  padding-left: 1.5rem; }
.custom-control-inline {
  display: inline-flex;
  margin-right: 1rem; }
.custom-control-input {
  position: absolute;
  z-index: -1;
  opacity: 0; }
  .custom-control-input:checked ~ .custom-control-label::before {
    color: #fff;
    background-color: #333; }
  .custom-control-input:focus ~ .custom-control-label::before {
    box-shadow: 0 0 0 1px #fff, 0 0 0 0.2rem rgba(51, 51, 51, 0.25); }
  .custom-control-input:active ~ .custom-control-label::before {
    color: #fff;
    background-color: #8c8c8c; }
  .custom-control-input:disabled ~ .custom-control-label {
    color: #6c757d; }
    .custom-control-input:disabled ~ .custom-control-label::before {
      background-color: #e9ecef; }
.custom-control-label {
  margin-bottom: 0; }
  .custom-control-label::before {
    position: absolute;
    top: 0.7rem;
    left: 0;
    display: block;
    width: 1rem;
    height: 1rem;
    pointer-events: none;
    content: "";
    user-select: none;
    background-color: #dee2e6; }
  .custom-control-label::after {
    position: absolute;
    top: 0.7rem;
    left: 0;
    display: block;
    width: 1rem;
    height: 1rem;
    content: "";
    background-repeat: no-repeat;
    background-position: center center;
    background-size: 50% 50%; }
.custom-checkbox .custom-control-label::before {
  border-radius: 0px; }
.custom-checkbox .custom-control-input:checked ~ .custom-control-label::before {
  background-color: #333; }
.custom-checkbox .custom-control-input:checked ~ .custom-control-label::after {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 8 8'%3E%3Cpath fill='%23fff' d='M6.564.75l-3.59 3.612-1.538-1.55L0 4.26 2.974 7.25 8 2.193z'/%3E%3C/svg%3E"); }
.custom-checkbox .custom-control-input:indeterminate ~ .custom-control-label::before {
  background-color: #333; }
.custom-checkbox .custom-control-input:indeterminate ~ .custom-control-label::after {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 4 4'%3E%3Cpath stroke='%23fff' d='M0 2h4'/%3E%3C/svg%3E"); }
.custom-checkbox .custom-control-input:disabled:checked ~ .custom-control-label::before {
  background-color: rgba(130, 76, 113, 0.5); }
.custom-checkbox .custom-control-input:disabled:indeterminate ~ .custom-control-label::before {
  background-color: rgba(130, 76, 113, 0.5); }
.custom-radio .custom-control-label::before {
  border-radius: 50%; }
.custom-radio .custom-control-input:checked ~ .custom-control-label::before {
  background-color: #333; }
.custom-radio .custom-control-input:checked ~ .custom-control-label::after {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='-4 -4 8 8'%3E%3Ccircle r='3' fill='%23fff'/%3E%3C/svg%3E"); }
.custom-radio .custom-control-input:disabled:checked ~ .custom-control-label::before {
  background-color: rgba(130, 76, 113, 0.5); }
.custom-select {
  display: inline-block;
  width: 100%;
  height: calc(3.15rem + 2px);
  padding: 0.375rem 1.75rem 0.375rem 0.75rem;
  line-height: 2.4;
  color: #495057;
  vertical-align: middle;
  background: #fff url("data:image/svg+xml;charset=utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 4 5'%3E%3Cpath fill='%23343a40' d='M2 0L0 2h4zm0 5L0 3h4z'/%3E%3C/svg%3E") no-repeat right 0.75rem center;
  background-size: 8px 10px;
  border: 1px solid #ced4da;
  border-radius: 0px;
  appearance: none; }
  .custom-select:focus {
    border-color: #737373;
    outline: 0;
    box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.075), 0 0 5px rgba(115, 115, 115, 0.5); }
    .custom-select:focus::-ms-value {
      color: #495057;
      background-color: #fff; }
  .custom-select[multiple], .custom-select[size]:not([size="1"]) {
    height: auto;
    padding-right: 0.75rem;
    background-image: none; }
  .custom-select:disabled {
    color: #6c757d;
    background-color: #e9ecef; }
  .custom-select::-ms-expand {
    opacity: 0; }
.custom-select-sm {
  height: calc(1.8125rem + 2px);
  padding-top: 0.375rem;
  padding-bottom: 0.375rem;
  font-size: 75%; }
.custom-select-lg {
  height: calc(2.875rem + 2px);
  padding-top: 0.375rem;
  padding-bottom: 0.375rem;
  font-size: 125%; }
.custom-file {
  position: relative;
  display: inline-block;
  width: 100%;
  height: calc(3.15rem + 2px);
  margin-bottom: 0; }
.custom-file-input {
  position: relative;
  z-index: 2;
  width: 100%;
  height: calc(3.15rem + 2px);
  margin: 0;
  opacity: 0; }
  .custom-file-input:focus ~ .custom-file-control {
    border-color: #737373;
    box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.25); }
    .custom-file-input:focus ~ .custom-file-control::before {
      border-color: #737373; }
  .custom-file-input:lang(en) ~ .custom-file-label::after {
    content: "Browse"; }
.custom-file-label {
  position: absolute;
  top: 0;
  right: 0;
  left: 0;
  z-index: 1;
  height: calc(3.15rem + 2px);
  padding: 0.375rem 0.75rem;
  line-height: 2.4;
  color: #495057;
  background-color: #fff;
  border: 1px solid #ced4da;
  border-radius: 0px; }
  .custom-file-label::after {
    position: absolute;
    top: 0;
    right: 0;
    bottom: 0;
    z-index: 3;
    display: block;
    height: calc(calc(3.15rem + 2px) - 1px * 2);
    padding: 0.375rem 0.75rem;
    line-height: 2.4;
    color: #495057;
    content: "Browse";
    background-color: #e9ecef;
    border-left: 1px solid #ced4da;
    border-radius: 0 0px 0px 0; }
.nav {
  display: flex;
  flex-wrap: wrap;
  padding-left: 0;
  margin-bottom: 0;
  list-style: none; }
.nav-link {
  display: block;
  padding: 0.5rem 1rem; }
  .nav-link:hover, .nav-link:focus {
    text-decoration: none; }
  .nav-link.disabled {
    color: #6c757d; }
.nav-tabs {
  border-bottom: 1px solid #dee2e6; }
  .nav-tabs .nav-item {
    margin-bottom: -1px; }
  .nav-tabs .nav-link {
    border: 1px solid transparent;
    border-top-left-radius: 0px;
    border-top-right-radius: 0px; }
    .nav-tabs .nav-link:hover, .nav-tabs .nav-link:focus {
      border-color: #e9ecef #e9ecef #dee2e6; }
    .nav-tabs .nav-link.disabled {
      color: #6c757d;
      background-color: transparent;
      border-color: transparent; }
  .nav-tabs .nav-link.active,
  .nav-tabs .nav-item.show .nav-link {
    color: #495057;
    background-color: #fff;
    border-color: #dee2e6 #dee2e6 #fff; }
  .nav-tabs .dropdown-menu {
    margin-top: -1px;
    border-top-left-radius: 0;
    border-top-right-radius: 0; }
.nav-pills .nav-link {
  border-radius: 0px; }
.nav-pills .nav-link.active,
.nav-pills .show > .nav-link {
  color: #fff;
  background-color: #333; }
.nav-fill .nav-item {
  flex: 1 1 auto;
  text-align: center; }
.nav-justified .nav-item {
  flex-basis: 0;
  flex-grow: 1;
  text-align: center; }
.tab-content > .tab-pane {
  display: none; }
.tab-content > .active {
  display: block; }
.navbar {
  position: relative;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  justify-content: space-between;
  padding: 0.5rem 1rem; }
  .navbar > .container,
  .navbar > .container-fluid {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between; }
.navbar-brand {
  display: inline-block;
  padding-top: 0.2rem;
  padding-bottom: 0.2rem;
  margin-right: 1rem;
  font-size: 1.25rem;
  line-height: inherit;
  white-space: nowrap; }
  .navbar-brand:hover, .navbar-brand:focus {
    text-decoration: none; }
.navbar-nav {
  display: flex;
  flex-direction: column;
  padding-left: 0;
  margin-bottom: 0;
  list-style: none; }
  .navbar-nav .nav-link {
    padding-right: 0;
    padding-left: 0; }
  .navbar-nav .dropdown-menu {
    position: static;
    float: none; }
.navbar-text {
  display: inline-block;
  padding-top: 0.5rem;
  padding-bottom: 0.5rem; }
.navbar-collapse {
  flex-basis: 100%;
  flex-grow: 1;
  align-items: center; }
.navbar-toggler {
  padding: 0.25rem 0.75rem;
  font-size: 1.25rem;
  line-height: 1;
  background-color: transparent;
  border: 1px solid transparent;
  border-radius: 0px; }
  .navbar-toggler:hover, .navbar-toggler:focus {
    text-decoration: none; }
  .navbar-toggler:not(:disabled):not(.disabled) {
    cursor: pointer; }
.navbar-toggler-icon {
  display: inline-block;
  width: 1.5em;
  height: 1.5em;
  vertical-align: middle;
  content: "";
  background: no-repeat center center;
  background-size: 100% 100%; }
@media (max-width: 575.98px) {
  .navbar-expand-sm > .container,
  .navbar-expand-sm > .container-fluid {
    padding-right: 0;
    padding-left: 0; } }
@media (min-width: 576px) {
  .navbar-expand-sm {
    flex-flow: row nowrap;
    justify-content: flex-start; }
    .navbar-expand-sm .navbar-nav {
      flex-direction: row; }
      .navbar-expand-sm .navbar-nav .dropdown-menu {
        position: absolute; }
      .navbar-expand-sm .navbar-nav .dropdown-menu-right {
        right: 0;
        left: auto; }
      .navbar-expand-sm .navbar-nav .nav-link {
        padding-right: 0.5rem;
        padding-left: 0.5rem; }
    .navbar-expand-sm > .container,
    .navbar-expand-sm > .container-fluid {
      flex-wrap: nowrap; }
    .navbar-expand-sm .navbar-collapse {
      display: flex !important;
      flex-basis: auto; }
    .navbar-expand-sm .navbar-toggler {
      display: none; }
    .navbar-expand-sm .dropup .dropdown-menu {
      top: auto;
      bottom: 100%; } }
@media (max-width: 767.98px) {
  .navbar-expand-md > .container,
  .navbar-expand-md > .container-fluid {
    padding-right: 0;
    padding-left: 0; } }
@media (min-width: 768px) {
  .navbar-expand-md {
    flex-flow: row nowrap;
    justify-content: flex-start; }
    .navbar-expand-md .navbar-nav {
      flex-direction: row; }
      .navbar-expand-md .navbar-nav .dropdown-menu {
        position: absolute; }
      .navbar-expand-md .navbar-nav .dropdown-menu-right {
        right: 0;
        left: auto; }
      .navbar-expand-md .navbar-nav .nav-link {
        padding-right: 0.5rem;
        padding-left: 0.5rem; }
    .navbar-expand-md > .container,
    .navbar-expand-md > .container-fluid {
      flex-wrap: nowrap; }
    .navbar-expand-md .navbar-collapse {
      display: flex !important;
      flex-basis: auto; }
    .navbar-expand-md .navbar-toggler {
      display: none; }
    .navbar-expand-md .dropup .dropdown-menu {
      top: auto;
      bottom: 100%; } }
@media (max-width: 991.98px) {
  .navbar-expand-lg > .container,
  .navbar-expand-lg > .container-fluid {
    padding-right: 0;
    padding-left: 0; } }
@media (min-width: 992px) {
  .navbar-expand-lg {
    flex-flow: row nowrap;
    justify-content: flex-start; }
    .navbar-expand-lg .navbar-nav {
      flex-direction: row; }
      .navbar-expand-lg .navbar-nav .dropdown-menu {
        position: absolute; }
      .navbar-expand-lg .navbar-nav .dropdown-menu-right {
        right: 0;
        left: auto; }
      .navbar-expand-lg .navbar-nav .nav-link {
        padding-right: 0.5rem;
        padding-left: 0.5rem; }
    .navbar-expand-lg > .container,
    .navbar-expand-lg > .container-fluid {
      flex-wrap: nowrap; }
    .navbar-expand-lg .navbar-collapse {
      display: flex !important;
      flex-basis: auto; }
    .navbar-expand-lg .navbar-toggler {
      display: none; }
    .navbar-expand-lg .dropup .dropdown-menu {
      top: auto;
      bottom: 100%; } }
@media (max-width: 1199.98px) {
  .navbar-expand-xl > .container,
  .navbar-expand-xl > .container-fluid {
    padding-right: 0;
    padding-left: 0; } }
@media (min-width: 1200px) {
  .navbar-expand-xl {
    flex-flow: row nowrap;
    justify-content: flex-start; }
    .navbar-expand-xl .navbar-nav {
      flex-direction: row; }
      .navbar-expand-xl .navbar-nav .dropdown-menu {
        position: absolute; }
      .navbar-expand-xl .navbar-nav .dropdown-menu-right {
        right: 0;
        left: auto; }
      .navbar-expand-xl .navbar-nav .nav-link {
        padding-right: 0.5rem;
        padding-left: 0.5rem; }
    .navbar-expand-xl > .container,
    .navbar-expand-xl > .container-fluid {
      flex-wrap: nowrap; }
    .navbar-expand-xl .navbar-collapse {
      display: flex !important;
      flex-basis: auto; }
    .navbar-expand-xl .navbar-toggler {
      display: none; }
    .navbar-expand-xl .dropup .dropdown-menu {
      top: auto;
      bottom: 100%; } }
.navbar-expand {
  flex-flow: row nowrap;
  justify-content: flex-start; }
  .navbar-expand > .container,
  .navbar-expand > .container-fluid {
    padding-right: 0;
    padding-left: 0; }
  .navbar-expand .navbar-nav {
    flex-direction: row; }
    .navbar-expand .navbar-nav .dropdown-menu {
      position: absolute; }
    .navbar-expand .navbar-nav .dropdown-menu-right {
      right: 0;
      left: auto; }
    .navbar-expand .navbar-nav .nav-link {
      padding-right: 0.5rem;
      padding-left: 0.5rem; }
  .navbar-expand > .container,
  .navbar-expand > .container-fluid {
    flex-wrap: nowrap; }
  .navbar-expand .navbar-collapse {
    display: flex !important;
    flex-basis: auto; }
  .navbar-expand .navbar-toggler {
    display: none; }
  .navbar-expand .dropup .dropdown-menu {
    top: auto;
    bottom: 100%; }
.navbar-light .navbar-brand {
  color: rgba(0, 0, 0, 0.9); }
  .navbar-light .navbar-brand:hover, .navbar-light .navbar-brand:focus {
    color: rgba(0, 0, 0, 0.9); }
.navbar-light .navbar-nav .nav-link {
  color: rgba(0, 0, 0, 0.5); }
  .navbar-light .navbar-nav .nav-link:hover, .navbar-light .navbar-nav .nav-link:focus {
    color: rgba(0, 0, 0, 0.7); }
  .navbar-light .navbar-nav .nav-link.disabled {
    color: rgba(0, 0, 0, 0.3); }
.navbar-light .navbar-nav .show > .nav-link,
.navbar-light .navbar-nav .active > .nav-link,
.navbar-light .navbar-nav .nav-link.show,
.navbar-light .navbar-nav .nav-link.active {
  color: rgba(0, 0, 0, 0.9); }
.navbar-light .navbar-toggler {
  color: rgba(0, 0, 0, 0.5);
  border-color: rgba(0, 0, 0, 0.1); }
.navbar-light .navbar-toggler-icon {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg viewBox='0 0 30 30' xmlns='http://www.w3.org/2000/svg'%3E%3Cpath stroke='rgba(0, 0, 0, 0.5)' stroke-width='2' stroke-linecap='round' stroke-miterlimit='10' d='M4 7h22M4 15h22M4 23h22'/%3E%3C/svg%3E"); }
.navbar-light .navbar-text {
  color: rgba(0, 0, 0, 0.5); }
  .navbar-light .navbar-text a {
    color: rgba(0, 0, 0, 0.9); }
    .navbar-light .navbar-text a:hover, .navbar-light .navbar-text a:focus {
      color: rgba(0, 0, 0, 0.9); }
.navbar-dark .navbar-brand {
  color: #fff; }
  .navbar-dark .navbar-brand:hover, .navbar-dark .navbar-brand:focus {
    color: #fff; }
.navbar-dark .navbar-nav .nav-link {
  color: rgba(255, 255, 255, 0.5); }
  .navbar-dark .navbar-nav .nav-link:hover, .navbar-dark .navbar-nav .nav-link:focus {
    color: rgba(255, 255, 255, 0.75); }
  .navbar-dark .navbar-nav .nav-link.disabled {
    color: rgba(255, 255, 255, 0.25); }
.navbar-dark .navbar-nav .show > .nav-link,
.navbar-dark .navbar-nav .active > .nav-link,
.navbar-dark .navbar-nav .nav-link.show,
.navbar-dark .navbar-nav .nav-link.active {
  color: #fff; }
.navbar-dark .navbar-toggler {
  color: rgba(255, 255, 255, 0.5);
  border-color: rgba(255, 255, 255, 0.1); }
.navbar-dark .navbar-toggler-icon {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg viewBox='0 0 30 30' xmlns='http://www.w3.org/2000/svg'%3E%3Cpath stroke='rgba(255, 255, 255, 0.5)' stroke-width='2' stroke-linecap='round' stroke-miterlimit='10' d='M4 7h22M4 15h22M4 23h22'/%3E%3C/svg%3E"); }
.navbar-dark .navbar-text {
  color: rgba(255, 255, 255, 0.5); }
  .navbar-dark .navbar-text a {
    color: #fff; }
    .navbar-dark .navbar-text a:hover, .navbar-dark .navbar-text a:focus {
      color: #fff; }
.card {
  position: relative;
  display: flex;
  flex-direction: column;
  min-width: 0;
  word-wrap: break-word;
  background-color: #fff;
  background-clip: border-box;
  border: 1px solid rgba(0, 0, 0, 0.125);
  border-radius: 0px; }
  .card > hr {
    margin-right: 0;
    margin-left: 0; }
  .card > .list-group:first-child .list-group-item:first-child {
    border-top-left-radius: 0px;
    border-top-right-radius: 0px; }
  .card > .list-group:last-child .list-group-item:last-child {
    border-bottom-right-radius: 0px;
    border-bottom-left-radius: 0px; }
.card-body {
  flex: 1 1 auto;
  padding: 1.25rem; }
.card-title {
  margin-bottom: 0.75rem; }
.card-subtitle {
  margin-top: -0.375rem;
  margin-bottom: 0; }
.card-text:last-child {
  margin-bottom: 0; }
.card-link:hover {
  text-decoration: none; }
.card-link + .card-link {
  margin-left: 1.25rem; }
.card-header {
  padding: 0.75rem 1.25rem;
  margin-bottom: 0;
  background-color: rgba(0, 0, 0, 0.03);
  border-bottom: 1px solid rgba(0, 0, 0, 0.125); }
  .card-header:first-child {
    border-radius: calc(0px - 1px) calc(0px - 1px) 0 0; }
  .card-header + .list-group .list-group-item:first-child {
    border-top: 0; }
.card-footer {
  padding: 0.75rem 1.25rem;
  background-color: rgba(0, 0, 0, 0.03);
  border-top: 1px solid rgba(0, 0, 0, 0.125); }
  .card-footer:last-child {
    border-radius: 0 0 calc(0px - 1px) calc(0px - 1px); }
.card-header-tabs {
  margin-right: -0.625rem;
  margin-bottom: -0.75rem;
  margin-left: -0.625rem;
  border-bottom: 0; }
.card-header-pills {
  margin-right: -0.625rem;
  margin-left: -0.625rem; }
.card-img-overlay {
  position: absolute;
  top: 0;
  right: 0;
  bottom: 0;
  left: 0;
  padding: 1.25rem; }
.card-img {
  width: 100%;
  border-radius: calc(0px - 1px); }
.card-img-top {
  width: 100%;
  border-top-left-radius: calc(0px - 1px);
  border-top-right-radius: calc(0px - 1px); }
.card-img-bottom {
  width: 100%;
  border-bottom-right-radius: calc(0px - 1px);
  border-bottom-left-radius: calc(0px - 1px); }
.card-deck {
  display: flex;
  flex-direction: column; }
  .card-deck .card {
    margin-bottom: 15px; }
  @media (min-width: 576px) {
    .card-deck {
      flex-flow: row wrap;
      margin-right: -15px;
      margin-left: -15px; }
      .card-deck .card {
        display: flex;
        flex: 1 0 0%;
        flex-direction: column;
        margin-right: 15px;
        margin-bottom: 0;
        margin-left: 15px; } }
.card-group {
  display: flex;
  flex-direction: column; }
  .card-group > .card {
    margin-bottom: 15px; }
  @media (min-width: 576px) {
    .card-group {
      flex-flow: row wrap; }
      .card-group > .card {
        flex: 1 0 0%;
        margin-bottom: 0; }
        .card-group > .card + .card {
          margin-left: 0;
          border-left: 0; }
        .card-group > .card:first-child {
          border-top-right-radius: 0;
          border-bottom-right-radius: 0; }
          .card-group > .card:first-child .card-img-top,
          .card-group > .card:first-child .card-header {
            border-top-right-radius: 0; }
          .card-group > .card:first-child .card-img-bottom,
          .card-group > .card:first-child .card-footer {
            border-bottom-right-radius: 0; }
        .card-group > .card:last-child {
          border-top-left-radius: 0;
          border-bottom-left-radius: 0; }
          .card-group > .card:last-child .card-img-top,
          .card-group > .card:last-child .card-header {
            border-top-left-radius: 0; }
          .card-group > .card:last-child .card-img-bottom,
          .card-group > .card:last-child .card-footer {
            border-bottom-left-radius: 0; }
        .card-group > .card:only-child {
          border-radius: 0px; }
          .card-group > .card:only-child .card-img-top,
          .card-group > .card:only-child .card-header {
            border-top-left-radius: 0px;
            border-top-right-radius: 0px; }
          .card-group > .card:only-child .card-img-bottom,
          .card-group > .card:only-child .card-footer {
            border-bottom-right-radius: 0px;
            border-bottom-left-radius: 0px; }
        .card-group > .card:not(:first-child):not(:last-child):not(:only-child) {
          border-radius: 0; }
          .card-group > .card:not(:first-child):not(:last-child):not(:only-child) .card-img-top,
          .card-group > .card:not(:first-child):not(:last-child):not(:only-child) .card-img-bottom,
          .card-group > .card:not(:first-child):not(:last-child):not(:only-child) .card-header,
          .card-group > .card:not(:first-child):not(:last-child):not(:only-child) .card-footer {
            border-radius: 0; } }
.card-columns .card {
  margin-bottom: 0.75rem; }
@media (min-width: 576px) {
  .card-columns {
    column-count: 3;
    column-gap: 1.25rem; }
    .card-columns .card {
      display: inline-block;
      width: 100%; } }
.breadcrumb {
  display: flex;
  flex-wrap: wrap;
  padding: 0.75rem 1rem;
  margin-bottom: 1rem;
  list-style: none;
  background-color: #e9ecef;
  border-radius: 0px; }
.breadcrumb-item + .breadcrumb-item::before {
  display: inline-block;
  padding-right: 0.5rem;
  padding-left: 0.5rem;
  color: #6c757d;
  content: "/"; }
.breadcrumb-item + .breadcrumb-item:hover::before {
  text-decoration: underline; }
.breadcrumb-item + .breadcrumb-item:hover::before {
  text-decoration: none; }
.breadcrumb-item.active {
  color: #6c757d; }
.pagination {
  display: flex;
  padding-left: 0;
  list-style: none;
  border-radius: 0px; }
.page-link {
  position: relative;
  display: block;
  padding: 0.5rem 0.75rem;
  margin-left: -1px;
  line-height: 1.25;
  color: #512479;
  background-color: #fff;
  border: 1px solid #dee2e6; }
  .page-link:hover {
    color: #523047;
    text-decoration: none;
    background-color: #e9ecef;
    border-color: #dee2e6; }
  .page-link:focus {
    z-index: 2;
    outline: 0;
    box-shadow: 0 0 0 0.2rem rgba(51, 51, 51, 0.25); }
  .page-link:not(:disabled):not(.disabled) {
    cursor: pointer; }
.page-item:first-child .page-link {
  margin-left: 0;
  border-top-left-radius: 0px;
  border-bottom-left-radius: 0px; }
.page-item:last-child .page-link {
  border-top-right-radius: 0px;
  border-bottom-right-radius: 0px; }
.page-item.active .page-link {
  z-index: 1;
  color: #fff;
  background-color: #333;
  border-color: #333; }
.page-item.disabled .page-link {
  color: #6c757d;
  pointer-events: none;
  cursor: auto;
  background-color: #fff;
  border-color: #dee2e6; }
.pagination-lg .page-link {
  padding: 0.75rem 1.5rem;
  font-size: 1.25rem;
  line-height: 1.5; }
.pagination-lg .page-item:first-child .page-link {
  border-top-left-radius: 0px;
  border-bottom-left-radius: 0px; }
.pagination-lg .page-item:last-child .page-link {
  border-top-right-radius: 0px;
  border-bottom-right-radius: 0px; }
.pagination-sm .page-link {
  padding: 0.25rem 0.5rem;
  font-size: 0.875rem;
  line-height: 1.5; }
.pagination-sm .page-item:first-child .page-link {
  border-top-left-radius: 0px;
  border-bottom-left-radius: 0px; }
.pagination-sm .page-item:last-child .page-link {
  border-top-right-radius: 0px;
  border-bottom-right-radius: 0px; }
.badge {
  display: inline-block;
  padding: 0.25em 0.4em;
  font-size: 75%;
  font-weight: 700;
  line-height: 1;
  text-align: center;
  white-space: nowrap;
  vertical-align: baseline;
  border-radius: 0px; }
  .badge:empty {
    display: none; }
.btn .badge {
  position: relative;
  top: -1px; }
.badge-pill {
  padding-right: 0.6em;
  padding-left: 0.6em;
  border-radius: 10rem; }
.badge-primary {
  color: #fff;
  background-color: #512479; }
  .badge-primary[href]:hover, .badge-primary[href]:focus {
    color: #fff;
    text-decoration: none;
    background-color: #49075e; }
.badge-secondary {
  color: #fff;
  background-color: #333; }
  .badge-secondary[href]:hover, .badge-secondary[href]:focus {
    color: #fff;
    text-decoration: none;
    background-color: #1a1a1a; }
.badge-success {
  color: #fff;
  background-color: #28a745; }
  .badge-success[href]:hover, .badge-success[href]:focus {
    color: #fff;
    text-decoration: none;
    background-color: #1e7e34; }
.badge-info {
  color: #fff;
  background-color: #17a2b8; }
  .badge-info[href]:hover, .badge-info[href]:focus {
    color: #fff;
    text-decoration: none;
    background-color: #117a8b; }
.badge-warning {
  color: #212529;
  background-color: #ffc107; }
  .badge-warning[href]:hover, .badge-warning[href]:focus {
    color: #212529;
    text-decoration: none;
    background-color: #d39e00; }
.badge-danger {
  color: #fff;
  background-color: #dc3545; }
  .badge-danger[href]:hover, .badge-danger[href]:focus {
    color: #fff;
    text-decoration: none;
    background-color: #bd2130; }
.badge-light {
  color: #212529;
  background-color: #f8f9fa; }
  .badge-light[href]:hover, .badge-light[href]:focus {
    color: #212529;
    text-decoration: none;
    background-color: #dae0e5; }
.badge-dark {
  color: #fff;
  background-color: #343a40; }
  .badge-dark[href]:hover, .badge-dark[href]:focus {
    color: #fff;
    text-decoration: none;
    background-color: #1d2124; }
.jumbotron {
  padding: 2rem 1rem;
  margin-bottom: 2rem;
  background-color: #e9ecef;
  border-radius: 0px; }
  @media (min-width: 576px) {
    .jumbotron {
      padding: 4rem 2rem; } }
.jumbotron-fluid {
  padding-right: 0;
  padding-left: 0;
  border-radius: 0; }
.alert {
  position: relative;
  padding: 0.75rem 1.25rem;
  margin-bottom: 1rem;
  border: 1px solid transparent;
  border-radius: 0px; }
.alert-heading {
  color: inherit; }
.alert-link {
  font-weight: 700; }
.alert-dismissible {
  padding-right: 4rem; }
  .alert-dismissible .close {
    position: absolute;
    top: 0;
    right: 0;
    padding: 0.75rem 1.25rem;
    color: inherit; }
.alert-primary {
  color: #44283b;
  background-color: #e6dbe3;
  border-color: #dccdd7; }
  .alert-primary hr {
    border-top-color: #d2becb; }
  .alert-primary .alert-link {
    color: #24151f; }
.alert-secondary {
  color: #1b1b1b;
  background-color: #d6d6d6;
  border-color: #c6c6c6; }
  .alert-secondary hr {
    border-top-color: #b9b9b9; }
  .alert-secondary .alert-link {
    color: #020202; }
.alert-success {
  color: #155724;
  background-color: #d4edda;
  border-color: #c3e6cb; }
  .alert-success hr {
    border-top-color: #b1dfbb; }
  .alert-success .alert-link {
    color: #0b2e13; }
.alert-info {
  color: #0c5460;
  background-color: #d1ecf1;
  border-color: #bee5eb; }
  .alert-info hr {
    border-top-color: #abdde5; }
  .alert-info .alert-link {
    color: #062c33; }
.alert-warning {
  color: #856404;
  background-color: #fff3cd;
  border-color: #ffeeba; }
  .alert-warning hr {
    border-top-color: #ffe8a1; }
  .alert-warning .alert-link {
    color: #533f03; }
.alert-danger {
  color: #721c24;
  background-color: #f8d7da;
  border-color: #f5c6cb; }
  .alert-danger hr {
    border-top-color: #f1b0b7; }
  .alert-danger .alert-link {
    color: #491217; }
.alert-light {
  color: #818182;
  background-color: #fefefe;
  border-color: #fdfdfe; }
  .alert-light hr {
    border-top-color: #ececf6; }
  .alert-light .alert-link {
    color: #686868; }
.alert-dark {
  color: #1b1e21;
  background-color: #d6d8d9;
  border-color: #c6c8ca; }
  .alert-dark hr {
    border-top-color: #b9bbbe; }
  .alert-dark .alert-link {
    color: #040505; }
@keyframes progress-bar-stripes {
  from {
    background-position: 1rem 0; }
  to {
    background-position: 0 0; } }
.progress {
  display: flex;
  height: 1rem;
  overflow: hidden;
  font-size: 0.75rem;
  background-color: #e9ecef;
  border-radius: 0px; }
.progress-bar {
  display: flex;
  flex-direction: column;
  justify-content: center;
  color: #fff;
  text-align: center;
  background-color: #512479;
  transition: width 0.6s ease; }
.progress-bar-striped {
  background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.15) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.15) 50%, rgba(255, 255, 255, 0.15) 75%, transparent 75%, transparent);
  background-size: 1rem 1rem; }
.progress-bar-animated {
  animation: progress-bar-stripes 1s linear infinite; }
.media {
  display: flex;
  align-items: flex-start; }
.media-body {
  flex: 1; }
.list-group {
  display: flex;
  flex-direction: column;
  padding-left: 0;
  margin-bottom: 0; }
.list-group-item-action {
  width: 100%;
  color: #495057;
  text-align: inherit; }
  .list-group-item-action:hover, .list-group-item-action:focus {
    color: #495057;
    text-decoration: none;
    background-color: #f8f9fa; }
  .list-group-item-action:active {
    color: #999;
    background-color: #e9ecef; }
.list-group-item {
  position: relative;
  display: block;
  padding: 0.75rem 1.25rem;
  margin-bottom: -1px;
  background-color: #fff;
  border: 1px solid rgba(0, 0, 0, 0.125); }
  .list-group-item:first-child {
    border-top-left-radius: 0px;
    border-top-right-radius: 0px; }
  .list-group-item:last-child {
    margin-bottom: 0;
    border-bottom-right-radius: 0px;
    border-bottom-left-radius: 0px; }
  .list-group-item:hover, .list-group-item:focus {
    z-index: 1;
    text-decoration: none; }
  .list-group-item.disabled, .list-group-item:disabled {
    color: #6c757d;
    background-color: #fff; }
  .list-group-item.active {
    z-index: 2;
    color: #fff;
    background-color: #333;
    border-color: #333; }
.list-group-flush .list-group-item {
  border-right: 0;
  border-left: 0;
  border-radius: 0; }
.list-group-flush:first-child .list-group-item:first-child {
  border-top: 0; }
.list-group-flush:last-child .list-group-item:last-child {
  border-bottom: 0; }
.list-group-item-primary {
  color: #44283b;
  background-color: #dccdd7; }
  .list-group-item-primary.list-group-item-action:hover, .list-group-item-primary.list-group-item-action:focus {
    color: #44283b;
    background-color: #d2becb; }
  .list-group-item-primary.list-group-item-action.active {
    color: #fff;
    background-color: #44283b;
    border-color: #44283b; }
.list-group-item-secondary {
  color: #1b1b1b;
  background-color: #c6c6c6; }
  .list-group-item-secondary.list-group-item-action:hover, .list-group-item-secondary.list-group-item-action:focus {
    color: #1b1b1b;
    background-color: #b9b9b9; }
  .list-group-item-secondary.list-group-item-action.active {
    color: #fff;
    background-color: #1b1b1b;
    border-color: #1b1b1b; }
.list-group-item-success {
  color: #155724;
  background-color: #c3e6cb; }
  .list-group-item-success.list-group-item-action:hover, .list-group-item-success.list-group-item-action:focus {
    color: #155724;
    background-color: #b1dfbb; }
  .list-group-item-success.list-group-item-action.active {
    color: #fff;
    background-color: #155724;
    border-color: #155724; }
.list-group-item-info {
  color: #0c5460;
  background-color: #bee5eb; }
  .list-group-item-info.list-group-item-action:hover, .list-group-item-info.list-group-item-action:focus {
    color: #0c5460;
    background-color: #abdde5; }
  .list-group-item-info.list-group-item-action.active {
    color: #fff;
    background-color: #0c5460;
    border-color: #0c5460; }
.list-group-item-warning {
  color: #856404;
  background-color: #ffeeba; }
  .list-group-item-warning.list-group-item-action:hover, .list-group-item-warning.list-group-item-action:focus {
    color: #856404;
    background-color: #ffe8a1; }
  .list-group-item-warning.list-group-item-action.active {
    color: #fff;
    background-color: #856404;
    border-color: #856404; }
.list-group-item-danger {
  color: #721c24;
  background-color: #f5c6cb; }
  .list-group-item-danger.list-group-item-action:hover, .list-group-item-danger.list-group-item-action:focus {
    color: #721c24;
    background-color: #f1b0b7; }
  .list-group-item-danger.list-group-item-action.active {
    color: #fff;
    background-color: #721c24;
    border-color: #721c24; }
.list-group-item-light {
  color: #818182;
  background-color: #fdfdfe; }
  .list-group-item-light.list-group-item-action:hover, .list-group-item-light.list-group-item-action:focus {
    color: #818182;
    background-color: #ececf6; }
  .list-group-item-light.list-group-item-action.active {
    color: #fff;
    background-color: #818182;
    border-color: #818182; }
.list-group-item-dark {
  color: #1b1e21;
  background-color: #c6c8ca; }
  .list-group-item-dark.list-group-item-action:hover, .list-group-item-dark.list-group-item-action:focus {
    color: #1b1e21;
    background-color: #b9bbbe; }
  .list-group-item-dark.list-group-item-action.active {
    color: #fff;
    background-color: #1b1e21;
    border-color: #1b1e21; }
.close {
  float: right;
  font-size: 1.5rem;
  font-weight: 700;
  line-height: 1;
  color: #000;
  text-shadow: 0 1px 0 #fff;
  opacity: .5; }
  .close:hover, .close:focus {
    color: #000;
    text-decoration: none;
    opacity: .75; }
  .close:not(:disabled):not(.disabled) {
    cursor: pointer; }
button.close {
  padding: 0;
  background-color: transparent;
  border: 0;
  -webkit-appearance: none; }
.modal-open {
  overflow: hidden; }
.modal {
  position: fixed;
  top: 0;
  right: 0;
  bottom: 0;
  left: 0;
  z-index: 1050;
  display: none;
  overflow: hidden;
  outline: 0; }
  .modal-open .modal {
    overflow-x: hidden;
    overflow-y: auto; }
.modal-dialog {
  position: relative;
  width: auto;
  margin: 0.5rem;
  pointer-events: none; }
  .modal.fade .modal-dialog {
    transition: transform 0.3s ease-out;
    transform: translate(0, -25%); }
  .modal.show .modal-dialog {
    transform: translate(0, 0); }
.modal-dialog-centered {
  display: flex;
  align-items: center;
  min-height: calc(100% - (0.5rem * 2)); }
.modal-content {
  position: relative;
  display: flex;
  flex-direction: column;
  width: 100%;
  pointer-events: auto;
  background-color: #fff;
  background-clip: padding-box;
  border: 1px solid rgba(0, 0, 0, 0.2);
  border-radius: 0px;
  outline: 0; }
.modal-backdrop {
  position: fixed;
  top: 0;
  right: 0;
  bottom: 0;
  left: 0;
  z-index: 1040;
  background-color: #000; }
  .modal-backdrop.fade {
    opacity: 0; }
  .modal-backdrop.show {
    opacity: 0.5; }
.modal-header {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  padding: 1rem;
  border-bottom: 1px solid #e9ecef;
  border-top-left-radius: 0px;
  border-top-right-radius: 0px; }
  .modal-header .close {
    padding: 1rem;
    margin: -1rem -1rem -1rem auto; }
.modal-title {
  margin-bottom: 0;
  line-height: 2.4; }
.modal-body {
  position: relative;
  flex: 1 1 auto;
  padding: 1rem; }
.modal-footer {
  display: flex;
  align-items: center;
  justify-content: flex-end;
  padding: 1rem;
  border-top: 1px solid #e9ecef; }
  .modal-footer > :not(:first-child) {
    margin-left: .25rem; }
  .modal-footer > :not(:last-child) {
    margin-right: .25rem; }
.modal-scrollbar-measure {
  position: absolute;
  top: -9999px;
  width: 50px;
  height: 50px;
  overflow: scroll; }
@media (min-width: 576px) {
  .modal-dialog {
    max-width: 500px;
    margin: 1.75rem auto; }
  .modal-dialog-centered {
    min-height: calc(100% - (1.75rem * 2)); }
  .modal-sm {
    max-width: 300px; } }
@media (min-width: 992px) {
  .modal-lg {
    max-width: 800px; } }
.tooltip {
  position: absolute;
  z-index: 1070;
  display: block;
  margin: 0;
  font-family: "Poppins", sans-serif;
  font-style: normal;
  font-weight: 400;
  line-height: 2.4;
  text-align: left;
  text-align: start;
  text-decoration: none;
  text-shadow: none;
  text-transform: none;
  letter-spacing: normal;
  word-break: normal;
  word-spacing: normal;
  white-space: normal;
  line-break: auto;
  font-size: 0.875rem;
  word-wrap: break-word;
  opacity: 0; }
  .tooltip.show {
    opacity: 0.9; }
  .tooltip .arrow {
    position: absolute;
    display: block;
    width: 0.8rem;
    height: 0.4rem; }
    .tooltip .arrow::before {
      position: absolute;
      content: "";
      border-color: transparent;
      border-style: solid; }
.bs-tooltip-top, .bs-tooltip-auto[x-placement^="top"] {
  padding: 0.4rem 0; }
  .bs-tooltip-top .arrow, .bs-tooltip-auto[x-placement^="top"] .arrow {
    bottom: 0; }
    .bs-tooltip-top .arrow::before, .bs-tooltip-auto[x-placement^="top"] .arrow::before {
      top: 0;
      border-width: 0.4rem 0.4rem 0;
      border-top-color: #000; }
.bs-tooltip-right, .bs-tooltip-auto[x-placement^="right"] {
  padding: 0 0.4rem; }
  .bs-tooltip-right .arrow, .bs-tooltip-auto[x-placement^="right"] .arrow {
    left: 0;
    width: 0.4rem;
    height: 0.8rem; }
    .bs-tooltip-right .arrow::before, .bs-tooltip-auto[x-placement^="right"] .arrow::before {
      right: 0;
      border-width: 0.4rem 0.4rem 0.4rem 0;
      border-right-color: #000; }
.bs-tooltip-bottom, .bs-tooltip-auto[x-placement^="bottom"] {
  padding: 0.4rem 0; }
  .bs-tooltip-bottom .arrow, .bs-tooltip-auto[x-placement^="bottom"] .arrow {
    top: 0; }
    .bs-tooltip-bottom .arrow::before, .bs-tooltip-auto[x-placement^="bottom"] .arrow::before {
      bottom: 0;
      border-width: 0 0.4rem 0.4rem;
      border-bottom-color: #000; }
.bs-tooltip-left, .bs-tooltip-auto[x-placement^="left"] {
  padding: 0 0.4rem; }
  .bs-tooltip-left .arrow, .bs-tooltip-auto[x-placement^="left"] .arrow {
    right: 0;
    width: 0.4rem;
    height: 0.8rem; }
    .bs-tooltip-left .arrow::before, .bs-tooltip-auto[x-placement^="left"] .arrow::before {
      left: 0;
      border-width: 0.4rem 0 0.4rem 0.4rem;
      border-left-color: #000; }
.tooltip-inner {
  max-width: 200px;
  padding: 0.25rem 0.5rem;
  color: #fff;
  text-align: center;
  background-color: #000;
  border-radius: 0px; }
.popover {
  position: absolute;
  top: 0;
  left: 0;
  z-index: 1060;
  display: block;
  max-width: 276px;
  font-family: "Poppins", sans-serif;
  font-style: normal;
  font-weight: 400;
  line-height: 2.4;
  text-align: left;
  text-align: start;
  text-decoration: none;
  text-shadow: none;
  text-transform: none;
  letter-spacing: normal;
  word-break: normal;
  word-spacing: normal;
  white-space: normal;
  line-break: auto;
  font-size: 0.875rem;
  word-wrap: break-word;
  background-color: #fff;
  background-clip: padding-box;
  border: 1px solid rgba(0, 0, 0, 0.2);
  border-radius: 0px; }
  .popover .arrow {
    position: absolute;
    display: block;
    width: 1rem;
    height: 0.5rem;
    margin: 0 0px; }
    .popover .arrow::before, .popover .arrow::after {
      position: absolute;
      display: block;
      content: "";
      border-color: transparent;
      border-style: solid; }
.bs-popover-top, .bs-popover-auto[x-placement^="top"] {
  margin-bottom: 0.5rem; }
  .bs-popover-top .arrow, .bs-popover-auto[x-placement^="top"] .arrow {
    bottom: calc((0.5rem + 1px) * -1); }
  .bs-popover-top .arrow::before, .bs-popover-auto[x-placement^="top"] .arrow::before,
  .bs-popover-top .arrow::after,
  .bs-popover-auto[x-placement^="top"] .arrow::after {
    border-width: 0.5rem 0.5rem 0; }
  .bs-popover-top .arrow::before, .bs-popover-auto[x-placement^="top"] .arrow::before {
    bottom: 0;
    border-top-color: rgba(0, 0, 0, 0.25); }
  .bs-popover-top .arrow::after, .bs-popover-auto[x-placement^="top"] .arrow::after {
    bottom: 1px;
    border-top-color: #fff; }
.bs-popover-right, .bs-popover-auto[x-placement^="right"] {
  margin-left: 0.5rem; }
  .bs-popover-right .arrow, .bs-popover-auto[x-placement^="right"] .arrow {
    left: calc((0.5rem + 1px) * -1);
    width: 0.5rem;
    height: 1rem;
    margin: 0px 0; }
  .bs-popover-right .arrow::before, .bs-popover-auto[x-placement^="right"] .arrow::before,
  .bs-popover-right .arrow::after,
  .bs-popover-auto[x-placement^="right"] .arrow::after {
    border-width: 0.5rem 0.5rem 0.5rem 0; }
  .bs-popover-right .arrow::before, .bs-popover-auto[x-placement^="right"] .arrow::before {
    left: 0;
    border-right-color: rgba(0, 0, 0, 0.25); }
  .bs-popover-right .arrow::after, .bs-popover-auto[x-placement^="right"] .arrow::after {
    left: 1px;
    border-right-color: #fff; }
.bs-popover-bottom, .bs-popover-auto[x-placement^="bottom"] {
  margin-top: 0.5rem; }
  .bs-popover-bottom .arrow, .bs-popover-auto[x-placement^="bottom"] .arrow {
    top: calc((0.5rem + 1px) * -1); }
  .bs-popover-bottom .arrow::before, .bs-popover-auto[x-placement^="bottom"] .arrow::before,
  .bs-popover-bottom .arrow::after,
  .bs-popover-auto[x-placement^="bottom"] .arrow::after {
    border-width: 0 0.5rem 0.5rem 0.5rem; }
  .bs-popover-bottom .arrow::before, .bs-popover-auto[x-placement^="bottom"] .arrow::before {
    top: 0;
    border-bottom-color: rgba(0, 0, 0, 0.25); }
  .bs-popover-bottom .arrow::after, .bs-popover-auto[x-placement^="bottom"] .arrow::after {
    top: 1px;
    border-bottom-color: #fff; }
  .bs-popover-bottom .popover-header::before, .bs-popover-auto[x-placement^="bottom"] .popover-header::before {
    position: absolute;
    top: 0;
    left: 50%;
    display: block;
    width: 1rem;
    margin-left: -0.5rem;
    content: "";
    border-bottom: 1px solid #f7f7f7; }
.bs-popover-left, .bs-popover-auto[x-placement^="left"] {
  margin-right: 0.5rem; }
  .bs-popover-left .arrow, .bs-popover-auto[x-placement^="left"] .arrow {
    right: calc((0.5rem + 1px) * -1);
    width: 0.5rem;
    height: 1rem;
    margin: 0px 0; }
  .bs-popover-left .arrow::before, .bs-popover-auto[x-placement^="left"] .arrow::before,
  .bs-popover-left .arrow::after,
  .bs-popover-auto[x-placement^="left"] .arrow::after {
    border-width: 0.5rem 0 0.5rem 0.5rem; }
  .bs-popover-left .arrow::before, .bs-popover-auto[x-placement^="left"] .arrow::before {
    right: 0;
    border-left-color: rgba(0, 0, 0, 0.25); }
  .bs-popover-left .arrow::after, .bs-popover-auto[x-placement^="left"] .arrow::after {
    right: 1px;
    border-left-color: #fff; }
.popover-header {
  padding: 0.5rem 0.75rem;
  margin-bottom: 0;
  font-size: 1rem;
  color: inherit;
  background-color: #f7f7f7;
  border-bottom: 1px solid #ebebeb;
  border-top-left-radius: calc(0px - 1px);
  border-top-right-radius: calc(0px - 1px); }
  .popover-header:empty {
    display: none; }
.popover-body {
  padding: 0.5rem 0.75rem;
  color: #999; }
.carousel {
  position: relative; }
.carousel-inner {
  position: relative;
  width: 100%;
  overflow: hidden; }
.carousel-item {
  position: relative;
  display: none;
  align-items: center;
  width: 100%;
  transition: transform 0.6s ease;
  backface-visibility: hidden;
  perspective: 1000px; }
.carousel-item.active,
.carousel-item-next,
.carousel-item-prev {
  display: block; }
.carousel-item-next,
.carousel-item-prev {
  position: absolute;
  top: 0; }
.carousel-item-next.carousel-item-left,
.carousel-item-prev.carousel-item-right {
  transform: translateX(0); }
  @supports (transform-style: preserve-3d) {
    .carousel-item-next.carousel-item-left,
    .carousel-item-prev.carousel-item-right {
      transform: translate3d(0, 0, 0); } }
.carousel-item-next,
.active.carousel-item-right {
  transform: translateX(100%); }
  @supports (transform-style: preserve-3d) {
    .carousel-item-next,
    .active.carousel-item-right {
      transform: translate3d(100%, 0, 0); } }
.carousel-item-prev,
.active.carousel-item-left {
  transform: translateX(-100%); }
  @supports (transform-style: preserve-3d) {
    .carousel-item-prev,
    .active.carousel-item-left {
      transform: translate3d(-100%, 0, 0); } }
.carousel-control-prev,
.carousel-control-next {
  position: absolute;
  top: 0;
  bottom: 0;
  display: flex;
  align-items: center;
  justify-content: center;
  width: 2%;
  color: #512479;
  text-align: center;
  opacity: 1; }
  .carousel-control-prev:hover, .carousel-control-prev:focus,
  .carousel-control-next:hover,
  .carousel-control-next:focus {
    color: #512479;
    text-decoration: none;
    outline: 0;
    opacity: .9; }
.carousel-control-prev {
  left: 0; }
.carousel-control-next {
  right: 0; }
.carousel-control-prev-icon,
.carousel-control-next-icon {
  display: inline-block;
  width: 20px;
  height: 20px;
  background: transparent no-repeat center center;
  background-size: 100% 100%; }
.carousel-control-prev-icon {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' fill='%23512479' viewBox='0 0 8 8'%3E%3Cpath d='M5.25 0l-4 4 4 4 1.5-1.5-2.5-2.5 2.5-2.5-1.5-1.5z'/%3E%3C/svg%3E"); }
.carousel-control-next-icon {
  background-image: url("data:image/svg+xml;charset=utf8,%3Csvg xmlns='http://www.w3.org/2000/svg' fill='%23512479' viewBox='0 0 8 8'%3E%3Cpath d='M2.75 0l-1.5 1.5 2.5 2.5-2.5 2.5 1.5 1.5 4-4-4-4z'/%3E%3C/svg%3E"); }
.carousel-indicators {
  position: absolute;
  right: 0;
  bottom: 10px;
  left: 0;
  z-index: 15;
  display: flex;
  justify-content: center;
  padding-left: 0;
  margin-right: 2%;
  margin-left: 2%;
  list-style: none; }
  .carousel-indicators li {
    position: relative;
    flex: 0 1 auto;
    width: 30px;
    height: 3px;
    margin-right: 3px;
    margin-left: 3px;
    text-indent: -999px;
    background-color: rgba(255, 255, 255, 0.5); }
    .carousel-indicators li::before {
      position: absolute;
      top: -10px;
      left: 0;
      display: inline-block;
      width: 100%;
      height: 10px;
      content: ""; }
    .carousel-indicators li::after {
      position: absolute;
      bottom: -10px;
      left: 0;
      display: inline-block;
      width: 100%;
      height: 10px;
      content: ""; }
  .carousel-indicators .active {
    background-color: #fff; }
.carousel-caption {
  position: absolute;
  right: 15%;
  bottom: 20px;
  left: 15%;
  z-index: 10;
  padding-top: 20px;
  padding-bottom: 20px;
  color: #fff;
  text-align: center; }
.align-baseline {
  vertical-align: baseline !important; }
.align-top {
  vertical-align: top !important; }
.align-middle {
  vertical-align: middle !important; }
.align-bottom {
  vertical-align: bottom !important; }
.align-text-bottom {
  vertical-align: text-bottom !important; }
.align-text-top {
  vertical-align: text-top !important; }
.bg-primary {
  background-color: #512479 !important; }
a.bg-primary:hover, a.bg-primary:focus,
button.bg-primary:hover,
button.bg-primary:focus {
  background-color: #49075e !important; }
.bg-secondary {
  background-color: #333 !important; }
a.bg-secondary:hover, a.bg-secondary:focus,
button.bg-secondary:hover,
button.bg-secondary:focus {
  background-color: #1a1a1a !important; }
.bg-success {
  background-color: #28a745 !important; }
a.bg-success:hover, a.bg-success:focus,
button.bg-success:hover,
button.bg-success:focus {
  background-color: #1e7e34 !important; }
.bg-info {
  background-color: #17a2b8 !important; }
a.bg-info:hover, a.bg-info:focus,
button.bg-info:hover,
button.bg-info:focus {
  background-color: #117a8b !important; }
.bg-warning {
  background-color: #ffc107 !important; }
a.bg-warning:hover, a.bg-warning:focus,
button.bg-warning:hover,
button.bg-warning:focus {
  background-color: #d39e00 !important; }
.bg-danger {
  background-color: #dc3545 !important; }
a.bg-danger:hover, a.bg-danger:focus,
button.bg-danger:hover,
button.bg-danger:focus {
  background-color: #bd2130 !important; }
.bg-light {
  background-color: #f8f9fa !important; }
a.bg-light:hover, a.bg-light:focus,
button.bg-light:hover,
button.bg-light:focus {
  background-color: #dae0e5 !important; }
.bg-dark {
  background-color: #343a40 !important; }
a.bg-dark:hover, a.bg-dark:focus,
button.bg-dark:hover,
button.bg-dark:focus {
  background-color: #1d2124 !important; }
.bg-white {
  background-color: #fff !important; }
.bg-transparent {
  background-color: transparent !important; }
.border {
  border: 1px solid #dee2e6 !important; }
.border-top {
  border-top: 1px solid #dee2e6 !important; }
.border-right {
  border-right: 1px solid #dee2e6 !important; }
.border-bottom {
  border-bottom: 1px solid #dee2e6 !important; }
.border-left {
  border-left: 1px solid #dee2e6 !important; }
.border-0 {
  border: 0 !important; }
.border-top-0 {
  border-top: 0 !important; }
.border-right-0 {
  border-right: 0 !important; }
.border-bottom-0 {
  border-bottom: 0 !important; }
.border-left-0 {
  border-left: 0 !important; }
.border-primary {
  border-color: #512479 !important; }
.border-secondary {
  border-color: #333 !important; }
.border-success {
  border-color: #28a745 !important; }
.border-info {
  border-color: #17a2b8 !important; }
.border-warning {
  border-color: #ffc107 !important; }
.border-danger {
  border-color: #dc3545 !important; }
.border-light {
  border-color: #f8f9fa !important; }
.border-dark {
  border-color: #343a40 !important; }
.border-white {
  border-color: #fff !important; }
.rounded {
  border-radius: 0px !important; }
.rounded-top {
  border-top-left-radius: 0px !important;
  border-top-right-radius: 0px !important; }
.rounded-right {
  border-top-right-radius: 0px !important;
  border-bottom-right-radius: 0px !important; }
.rounded-bottom {
  border-bottom-right-radius: 0px !important;
  border-bottom-left-radius: 0px !important; }
.rounded-left {
  border-top-left-radius: 0px !important;
  border-bottom-left-radius: 0px !important; }
.rounded-circle {
  border-radius: 50% !important; }
.rounded-0 {
  border-radius: 0 !important; }
.clearfix::after {
  display: block;
  clear: both;
  content: ""; }
.d-none {
  display: none !important; }
.d-inline {
  display: inline !important; }
.d-inline-block {
  display: inline-block !important; }
.d-block {
  display: block !important; }
.d-table {
  display: table !important; }
.d-table-row {
  display: table-row !important; }
.d-table-cell {
  display: table-cell !important; }
.d-flex {
  display: flex !important; }
.d-inline-flex {
  display: inline-flex !important; }
@media (min-width: 576px) {
  .d-sm-none {
    display: none !important; }
  .d-sm-inline {
    display: inline !important; }
  .d-sm-inline-block {
    display: inline-block !important; }
  .d-sm-block {
    display: block !important; }
  .d-sm-table {
    display: table !important; }
  .d-sm-table-row {
    display: table-row !important; }
  .d-sm-table-cell {
    display: table-cell !important; }
  .d-sm-flex {
    display: flex !important; }
  .d-sm-inline-flex {
    display: inline-flex !important; } }
@media (min-width: 768px) {
  .d-md-none {
    display: none !important; }
  .d-md-inline {
    display: inline !important; }
  .d-md-inline-block {
    display: inline-block !important; }
  .d-md-block {
    display: block !important; }
  .d-md-table {
    display: table !important; }
  .d-md-table-row {
    display: table-row !important; }
  .d-md-table-cell {
    display: table-cell !important; }
  .d-md-flex {
    display: flex !important; }
  .d-md-inline-flex {
    display: inline-flex !important; } }
@media (min-width: 992px) {
  .d-lg-none {
    display: none !important; }
  .d-lg-inline {
    display: inline !important; }
  .d-lg-inline-block {
    display: inline-block !important; }
  .d-lg-block {
    display: block !important; }
  .d-lg-table {
    display: table !important; }
  .d-lg-table-row {
    display: table-row !important; }
  .d-lg-table-cell {
    display: table-cell !important; }
  .d-lg-flex {
    display: flex !important; }
  .d-lg-inline-flex {
    display: inline-flex !important; } }
@media (min-width: 1200px) {
  .d-xl-none {
    display: none !important; }
  .d-xl-inline {
    display: inline !important; }
  .d-xl-inline-block {
    display: inline-block !important; }
  .d-xl-block {
    display: block !important; }
  .d-xl-table {
    display: table !important; }
  .d-xl-table-row {
    display: table-row !important; }
  .d-xl-table-cell {
    display: table-cell !important; }
  .d-xl-flex {
    display: flex !important; }
  .d-xl-inline-flex {
    display: inline-flex !important; } }
@media print {
  .d-print-none {
    display: none !important; }
  .d-print-inline {
    display: inline !important; }
  .d-print-inline-block {
    display: inline-block !important; }
  .d-print-block {
    display: block !important; }
  .d-print-table {
    display: table !important; }
  .d-print-table-row {
    display: table-row !important; }
  .d-print-table-cell {
    display: table-cell !important; }
  .d-print-flex {
    display: flex !important; }
  .d-print-inline-flex {
    display: inline-flex !important; } }
.embed-responsive {
  position: relative;
  display: block;
  width: 100%;
  padding: 0;
  overflow: hidden; }
  .embed-responsive::before {
    display: block;
    content: ""; }
  .embed-responsive .embed-responsive-item,
  .embed-responsive iframe,
  .embed-responsive embed,
  .embed-responsive object,
  .embed-responsive video {
    position: absolute;
    top: 0;
    bottom: 0;
    left: 0;
    width: 100%;
    height: 100%;
    border: 0; }
.embed-responsive-21by9::before {
  padding-top: 42.8571428571%; }
.embed-responsive-16by9::before {
  padding-top: 56.25%; }
.embed-responsive-4by3::before {
  padding-top: 75%; }
.embed-responsive-1by1::before {
  padding-top: 100%; }
.flex-row {
  flex-direction: row !important; }
.flex-column {
  flex-direction: column !important; }
.flex-row-reverse {
  flex-direction: row-reverse !important; }
.flex-column-reverse {
  flex-direction: column-reverse !important; }
.flex-wrap {
  flex-wrap: wrap !important; }
.flex-nowrap {
  flex-wrap: nowrap !important; }
.flex-wrap-reverse {
  flex-wrap: wrap-reverse !important; }
.justify-content-start {
  justify-content: flex-start !important; }
.justify-content-end {
  justify-content: flex-end !important; }
.justify-content-center {
  justify-content: center !important; }
.justify-content-between {
  justify-content: space-between !important; }
.justify-content-around {
  justify-content: space-around !important; }
.align-items-start {
  align-items: flex-start !important; }
.align-items-end {
  align-items: flex-end !important; }
.align-items-center {
  align-items: center !important; }
.align-items-baseline {
  align-items: baseline !important; }
.align-items-stretch {
  align-items: stretch !important; }
.align-content-start {
  align-content: flex-start !important; }
.align-content-end {
  align-content: flex-end !important; }
.align-content-center {
  align-content: center !important; }
.align-content-between {
  align-content: space-between !important; }
.align-content-around {
  align-content: space-around !important; }
.align-content-stretch {
  align-content: stretch !important; }
.align-self-auto {
  align-self: auto !important; }
.align-self-start {
  align-self: flex-start !important; }
.align-self-end {
  align-self: flex-end !important; }
.align-self-center {
  align-self: center !important; }
.align-self-baseline {
  align-self: baseline !important; }
.align-self-stretch {
  align-self: stretch !important; }
@media (min-width: 576px) {
  .flex-sm-row {
    flex-direction: row !important; }
  .flex-sm-column {
    flex-direction: column !important; }
  .flex-sm-row-reverse {
    flex-direction: row-reverse !important; }
  .flex-sm-column-reverse {
    flex-direction: column-reverse !important; }
  .flex-sm-wrap {
    flex-wrap: wrap !important; }
  .flex-sm-nowrap {
    flex-wrap: nowrap !important; }
  .flex-sm-wrap-reverse {
    flex-wrap: wrap-reverse !important; }
  .justify-content-sm-start {
    justify-content: flex-start !important; }
  .justify-content-sm-end {
    justify-content: flex-end !important; }
  .justify-content-sm-center {
    justify-content: center !important; }
  .justify-content-sm-between {
    justify-content: space-between !important; }
  .justify-content-sm-around {
    justify-content: space-around !important; }
  .align-items-sm-start {
    align-items: flex-start !important; }
  .align-items-sm-end {
    align-items: flex-end !important; }
  .align-items-sm-center {
    align-items: center !important; }
  .align-items-sm-baseline {
    align-items: baseline !important; }
  .align-items-sm-stretch {
    align-items: stretch !important; }
  .align-content-sm-start {
    align-content: flex-start !important; }
  .align-content-sm-end {
    align-content: flex-end !important; }
  .align-content-sm-center {
    align-content: center !important; }
  .align-content-sm-between {
    align-content: space-between !important; }
  .align-content-sm-around {
    align-content: space-around !important; }
  .align-content-sm-stretch {
    align-content: stretch !important; }
  .align-self-sm-auto {
    align-self: auto !important; }
  .align-self-sm-start {
    align-self: flex-start !important; }
  .align-self-sm-end {
    align-self: flex-end !important; }
  .align-self-sm-center {
    align-self: center !important; }
  .align-self-sm-baseline {
    align-self: baseline !important; }
  .align-self-sm-stretch {
    align-self: stretch !important; } }
@media (min-width: 768px) {
  .flex-md-row {
    flex-direction: row !important; }
  .flex-md-column {
    flex-direction: column !important; }
  .flex-md-row-reverse {
    flex-direction: row-reverse !important; }
  .flex-md-column-reverse {
    flex-direction: column-reverse !important; }
  .flex-md-wrap {
    flex-wrap: wrap !important; }
  .flex-md-nowrap {
    flex-wrap: nowrap !important; }
  .flex-md-wrap-reverse {
    flex-wrap: wrap-reverse !important; }
  .justify-content-md-start {
    justify-content: flex-start !important; }
  .justify-content-md-end {
    justify-content: flex-end !important; }
  .justify-content-md-center {
    justify-content: center !important; }
  .justify-content-md-between {
    justify-content: space-between !important; }
  .justify-content-md-around {
    justify-content: space-around !important; }
  .align-items-md-start {
    align-items: flex-start !important; }
  .align-items-md-end {
    align-items: flex-end !important; }
  .align-items-md-center {
    align-items: center !important; }
  .align-items-md-baseline {
    align-items: baseline !important; }
  .align-items-md-stretch {
    align-items: stretch !important; }
  .align-content-md-start {
    align-content: flex-start !important; }
  .align-content-md-end {
    align-content: flex-end !important; }
  .align-content-md-center {
    align-content: center !important; }
  .align-content-md-between {
    align-content: space-between !important; }
  .align-content-md-around {
    align-content: space-around !important; }
  .align-content-md-stretch {
    align-content: stretch !important; }
  .align-self-md-auto {
    align-self: auto !important; }
  .align-self-md-start {
    align-self: flex-start !important; }
  .align-self-md-end {
    align-self: flex-end !important; }
  .align-self-md-center {
    align-self: center !important; }
  .align-self-md-baseline {
    align-self: baseline !important; }
  .align-self-md-stretch {
    align-self: stretch !important; } }
@media (min-width: 992px) {
  .flex-lg-row {
    flex-direction: row !important; }
  .flex-lg-column {
    flex-direction: column !important; }
  .flex-lg-row-reverse {
    flex-direction: row-reverse !important; }
  .flex-lg-column-reverse {
    flex-direction: column-reverse !important; }
  .flex-lg-wrap {
    flex-wrap: wrap !important; }
  .flex-lg-nowrap {
    flex-wrap: nowrap !important; }
  .flex-lg-wrap-reverse {
    flex-wrap: wrap-reverse !important; }
  .justify-content-lg-start {
    justify-content: flex-start !important; }
  .justify-content-lg-end {
    justify-content: flex-end !important; }
  .justify-content-lg-center {
    justify-content: center !important; }
  .justify-content-lg-between {
    justify-content: space-between !important; }
  .justify-content-lg-around {
    justify-content: space-around !important; }
  .align-items-lg-start {
    align-items: flex-start !important; }
  .align-items-lg-end {
    align-items: flex-end !important; }
  .align-items-lg-center {
    align-items: center !important; }
  .align-items-lg-baseline {
    align-items: baseline !important; }
  .align-items-lg-stretch {
    align-items: stretch !important; }
  .align-content-lg-start {
    align-content: flex-start !important; }
  .align-content-lg-end {
    align-content: flex-end !important; }
  .align-content-lg-center {
    align-content: center !important; }
  .align-content-lg-between {
    align-content: space-between !important; }
  .align-content-lg-around {
    align-content: space-around !important; }
  .align-content-lg-stretch {
    align-content: stretch !important; }
  .align-self-lg-auto {
    align-self: auto !important; }
  .align-self-lg-start {
    align-self: flex-start !important; }
  .align-self-lg-end {
    align-self: flex-end !important; }
  .align-self-lg-center {
    align-self: center !important; }
  .align-self-lg-baseline {
    align-self: baseline !important; }
  .align-self-lg-stretch {
    align-self: stretch !important; } }
@media (min-width: 1200px) {
  .flex-xl-row {
    flex-direction: row !important; }
  .flex-xl-column {
    flex-direction: column !important; }
  .flex-xl-row-reverse {
    flex-direction: row-reverse !important; }
  .flex-xl-column-reverse {
    flex-direction: column-reverse !important; }
  .flex-xl-wrap {
    flex-wrap: wrap !important; }
  .flex-xl-nowrap {
    flex-wrap: nowrap !important; }
  .flex-xl-wrap-reverse {
    flex-wrap: wrap-reverse !important; }
  .justify-content-xl-start {
    justify-content: flex-start !important; }
  .justify-content-xl-end {
    justify-content: flex-end !important; }
  .justify-content-xl-center {
    justify-content: center !important; }
  .justify-content-xl-between {
    justify-content: space-between !important; }
  .justify-content-xl-around {
    justify-content: space-around !important; }
  .align-items-xl-start {
    align-items: flex-start !important; }
  .align-items-xl-end {
    align-items: flex-end !important; }
  .align-items-xl-center {
    align-items: center !important; }
  .align-items-xl-baseline {
    align-items: baseline !important; }
  .align-items-xl-stretch {
    align-items: stretch !important; }
  .align-content-xl-start {
    align-content: flex-start !important; }
  .align-content-xl-end {
    align-content: flex-end !important; }
  .align-content-xl-center {
    align-content: center !important; }
  .align-content-xl-between {
    align-content: space-between !important; }
  .align-content-xl-around {
    align-content: space-around !important; }
  .align-content-xl-stretch {
    align-content: stretch !important; }
  .align-self-xl-auto {
    align-self: auto !important; }
  .align-self-xl-start {
    align-self: flex-start !important; }
  .align-self-xl-end {
    align-self: flex-end !important; }
  .align-self-xl-center {
    align-self: center !important; }
  .align-self-xl-baseline {
    align-self: baseline !important; }
  .align-self-xl-stretch {
    align-self: stretch !important; } }
.float-left {
  float: left !important; }
.float-right {
  float: right !important; }
.float-none {
  float: none !important; }
@media (min-width: 576px) {
  .float-sm-left {
    float: left !important; }
  .float-sm-right {
    float: right !important; }
  .float-sm-none {
    float: none !important; } }
@media (min-width: 768px) {
  .float-md-left {
    float: left !important; }
  .float-md-right {
    float: right !important; }
  .float-md-none {
    float: none !important; } }
@media (min-width: 992px) {
  .float-lg-left {
    float: left !important; }
  .float-lg-right {
    float: right !important; }
  .float-lg-none {
    float: none !important; } }
@media (min-width: 1200px) {
  .float-xl-left {
    float: left !important; }
  .float-xl-right {
    float: right !important; }
  .float-xl-none {
    float: none !important; } }
.position-static {
  position: static !important; }
.position-relative {
  position: relative !important; }
.position-absolute {
  position: absolute !important; }
.position-fixed {
  position: fixed !important; }
.position-sticky {
  position: sticky !important; }
.fixed-top {
  position: fixed;
  top: 0;
  right: 0;
  left: 0;
  z-index: 1030; }
.fixed-bottom {
  position: fixed;
  right: 0;
  bottom: 0;
  left: 0;
  z-index: 1030; }
@supports (position: sticky) {
  .sticky-top {
    position: sticky;
    top: 0;
    z-index: 1020; } }
.sr-only {
  position: absolute;
  width: 1px;
  height: 1px;
  padding: 0;
  overflow: hidden;
  clip: rect(0, 0, 0, 0);
  white-space: nowrap;
  clip-path: inset(50%);
  border: 0; }
.sr-only-focusable:active, .sr-only-focusable:focus {
  position: static;
  width: auto;
  height: auto;
  overflow: visible;
  clip: auto;
  white-space: normal;
  clip-path: none; }
.w-25 {
  width: 25% !important; }
.w-50 {
  width: 50% !important; }
.w-75 {
  width: 75% !important; }
.w-100 {
  width: 100% !important; }
.h-25 {
  height: 25% !important; }
.h-50 {
  height: 50% !important; }
.h-75 {
  height: 75% !important; }
.h-100 {
  height: 100% !important; }
.mw-100 {
  max-width: 100% !important; }
.mh-100 {
  max-height: 100% !important; }
.m-0 {
  margin: 0 !important; }
.mt-0,
.my-0 {
  margin-top: 0 !important; }
.mr-0,
.mx-0 {
  margin-right: 0 !important; }
.mb-0,
.my-0 {
  margin-bottom: 0 !important; }
.ml-0,
.mx-0 {
  margin-left: 0 !important; }
.m-1 {
  margin: 0.25rem !important; }
.mt-1,
.my-1 {
  margin-top: 0.25rem !important; }
.mr-1,
.mx-1 {
  margin-right: 0.25rem !important; }
.mb-1,
.my-1 {
  margin-bottom: 0.25rem !important; }
.ml-1,
.mx-1 {
  margin-left: 0.25rem !important; }
.m-2 {
  margin: 0.5rem !important; }
.mt-2,
.my-2 {
  margin-top: 0.5rem !important; }
.mr-2,
.mx-2 {
  margin-right: 0.5rem !important; }
.mb-2,
.my-2 {
  margin-bottom: 0.5rem !important; }
.ml-2,
.mx-2 {
  margin-left: 0.5rem !important; }
.m-3 {
  margin: 1rem !important; }
.mt-3,
.my-3 {
  margin-top: 1rem !important; }
.mr-3,
.mx-3 {
  margin-right: 1rem !important; }
.mb-3,
.my-3 {
  margin-bottom: 1rem !important; }
.ml-3,
.mx-3 {
  margin-left: 1rem !important; }
.m-4 {
  margin: 1.5rem !important; }
.mt-4,
.my-4 {
  margin-top: 1.5rem !important; }
.mr-4,
.mx-4 {
  margin-right: 1.5rem !important; }
.mb-4,
.my-4 {
  margin-bottom: 1.5rem !important; }
.ml-4,
.mx-4 {
  margin-left: 1.5rem !important; }
.m-5 {
  margin: 3rem !important; }
.mt-5,
.my-5 {
  margin-top: 3rem !important; }
.mr-5,
.mx-5 {
  margin-right: 3rem !important; }
.mb-5,
.my-5 {
  margin-bottom: 3rem !important; }
.ml-5,
.mx-5 {
  margin-left: 3rem !important; }
.p-0 {
  padding: 0 !important; }
.pt-0,
.py-0 {
  padding-top: 0 !important; }
.pr-0,
.px-0 {
  padding-right: 0 !important; }
.pb-0,
.py-0 {
  padding-bottom: 0 !important; }
.pl-0,
.px-0 {
  padding-left: 0 !important; }
.p-1 {
  padding: 0.25rem !important; }
.pt-1,
.py-1 {
  padding-top: 0.25rem !important; }
.pr-1,
.px-1 {
  padding-right: 0.25rem !important; }
.pb-1,
.py-1 {
  padding-bottom: 0.25rem !important; }
.pl-1,
.px-1 {
  padding-left: 0.25rem !important; }
.p-2 {
  padding: 0.5rem !important; }
.pt-2,
.py-2 {
  padding-top: 0.5rem !important; }
.pr-2,
.px-2 {
  padding-right: 0.5rem !important; }
.pb-2,
.py-2 {
  padding-bottom: 0.5rem !important; }
.pl-2,
.px-2 {
  padding-left: 0.5rem !important; }
.p-3 {
  padding: 1rem !important; }
.pt-3,
.py-3 {
  padding-top: 1rem !important; }
.pr-3,
.px-3 {
  padding-right: 1rem !important; }
.pb-3,
.py-3 {
  padding-bottom: 1rem !important; }
.pl-3,
.px-3 {
  padding-left: 1rem !important; }
.p-4 {
  padding: 1.5rem !important; }
.pt-4,
.py-4 {
  padding-top: 1.5rem !important; }
.pr-4,
.px-4 {
  padding-right: 1.5rem !important; }
.pb-4,
.py-4 {
  padding-bottom: 1.5rem !important; }
.pl-4,
.px-4 {
  padding-left: 1.5rem !important; }
.p-5 {
  padding: 3rem !important; }
.pt-5,
.py-5 {
  padding-top: 3rem !important; }
.pr-5,
.px-5 {
  padding-right: 3rem !important; }
.pb-5,
.py-5 {
  padding-bottom: 3rem !important; }
.pl-5,
.px-5 {
  padding-left: 3rem !important; }
.m-auto {
  margin: auto !important; }
.mt-auto,
.my-auto {
  margin-top: auto !important; }
.mr-auto,
.mx-auto {
  margin-right: auto !important; }
.mb-auto,
.my-auto {
  margin-bottom: auto !important; }
.ml-auto,
.mx-auto {
  margin-left: auto !important; }
@media (min-width: 576px) {
  .m-sm-0 {
    margin: 0 !important; }
  .mt-sm-0,
  .my-sm-0 {
    margin-top: 0 !important; }
  .mr-sm-0,
  .mx-sm-0 {
    margin-right: 0 !important; }
  .mb-sm-0,
  .my-sm-0 {
    margin-bottom: 0 !important; }
  .ml-sm-0,
  .mx-sm-0 {
    margin-left: 0 !important; }
  .m-sm-1 {
    margin: 0.25rem !important; }
  .mt-sm-1,
  .my-sm-1 {
    margin-top: 0.25rem !important; }
  .mr-sm-1,
  .mx-sm-1 {
    margin-right: 0.25rem !important; }
  .mb-sm-1,
  .my-sm-1 {
    margin-bottom: 0.25rem !important; }
  .ml-sm-1,
  .mx-sm-1 {
    margin-left: 0.25rem !important; }
  .m-sm-2 {
    margin: 0.5rem !important; }
  .mt-sm-2,
  .my-sm-2 {
    margin-top: 0.5rem !important; }
  .mr-sm-2,
  .mx-sm-2 {
    margin-right: 0.5rem !important; }
  .mb-sm-2,
  .my-sm-2 {
    margin-bottom: 0.5rem !important; }
  .ml-sm-2,
  .mx-sm-2 {
    margin-left: 0.5rem !important; }
  .m-sm-3 {
    margin: 1rem !important; }
  .mt-sm-3,
  .my-sm-3 {
    margin-top: 1rem !important; }
  .mr-sm-3,
  .mx-sm-3 {
    margin-right: 1rem !important; }
  .mb-sm-3,
  .my-sm-3 {
    margin-bottom: 1rem !important; }
  .ml-sm-3,
  .mx-sm-3 {
    margin-left: 1rem !important; }
  .m-sm-4 {
    margin: 1.5rem !important; }
  .mt-sm-4,
  .my-sm-4 {
    margin-top: 1.5rem !important; }
  .mr-sm-4,
  .mx-sm-4 {
    margin-right: 1.5rem !important; }
  .mb-sm-4,
  .my-sm-4 {
    margin-bottom: 1.5rem !important; }
  .ml-sm-4,
  .mx-sm-4 {
    margin-left: 1.5rem !important; }
  .m-sm-5 {
    margin: 3rem !important; }
  .mt-sm-5,
  .my-sm-5 {
    margin-top: 3rem !important; }
  .mr-sm-5,
  .mx-sm-5 {
    margin-right: 3rem !important; }
  .mb-sm-5,
  .my-sm-5 {
    margin-bottom: 3rem !important; }
  .ml-sm-5,
  .mx-sm-5 {
    margin-left: 3rem !important; }
  .p-sm-0 {
    padding: 0 !important; }
  .pt-sm-0,
  .py-sm-0 {
    padding-top: 0 !important; }
  .pr-sm-0,
  .px-sm-0 {
    padding-right: 0 !important; }
  .pb-sm-0,
  .py-sm-0 {
    padding-bottom: 0 !important; }
  .pl-sm-0,
  .px-sm-0 {
    padding-left: 0 !important; }
  .p-sm-1 {
    padding: 0.25rem !important; }
  .pt-sm-1,
  .py-sm-1 {
    padding-top: 0.25rem !important; }
  .pr-sm-1,
  .px-sm-1 {
    padding-right: 0.25rem !important; }
  .pb-sm-1,
  .py-sm-1 {
    padding-bottom: 0.25rem !important; }
  .pl-sm-1,
  .px-sm-1 {
    padding-left: 0.25rem !important; }
  .p-sm-2 {
    padding: 0.5rem !important; }
  .pt-sm-2,
  .py-sm-2 {
    padding-top: 0.5rem !important; }
  .pr-sm-2,
  .px-sm-2 {
    padding-right: 0.5rem !important; }
  .pb-sm-2,
  .py-sm-2 {
    padding-bottom: 0.5rem !important; }
  .pl-sm-2,
  .px-sm-2 {
    padding-left: 0.5rem !important; }
  .p-sm-3 {
    padding: 1rem !important; }
  .pt-sm-3,
  .py-sm-3 {
    padding-top: 1rem !important; }
  .pr-sm-3,
  .px-sm-3 {
    padding-right: 1rem !important; }
  .pb-sm-3,
  .py-sm-3 {
    padding-bottom: 1rem !important; }
  .pl-sm-3,
  .px-sm-3 {
    padding-left: 1rem !important; }
  .p-sm-4 {
    padding: 1.5rem !important; }
  .pt-sm-4,
  .py-sm-4 {
    padding-top: 1.5rem !important; }
  .pr-sm-4,
  .px-sm-4 {
    padding-right: 1.5rem !important; }
  .pb-sm-4,
  .py-sm-4 {
    padding-bottom: 1.5rem !important; }
  .pl-sm-4,
  .px-sm-4 {
    padding-left: 1.5rem !important; }
  .p-sm-5 {
    padding: 3rem !important; }
  .pt-sm-5,
  .py-sm-5 {
    padding-top: 3rem !important; }
  .pr-sm-5,
  .px-sm-5 {
    padding-right: 3rem !important; }
  .pb-sm-5,
  .py-sm-5 {
    padding-bottom: 3rem !important; }
  .pl-sm-5,
  .px-sm-5 {
    padding-left: 3rem !important; }
  .m-sm-auto {
    margin: auto !important; }
  .mt-sm-auto,
  .my-sm-auto {
    margin-top: auto !important; }
  .mr-sm-auto,
  .mx-sm-auto {
    margin-right: auto !important; }
  .mb-sm-auto,
  .my-sm-auto {
    margin-bottom: auto !important; }
  .ml-sm-auto,
  .mx-sm-auto {
    margin-left: auto !important; } }
@media (min-width: 768px) {
  .m-md-0 {
    margin: 0 !important; }
  .mt-md-0,
  .my-md-0 {
    margin-top: 0 !important; }
  .mr-md-0,
  .mx-md-0 {
    margin-right: 0 !important; }
  .mb-md-0,
  .my-md-0 {
    margin-bottom: 0 !important; }
  .ml-md-0,
  .mx-md-0 {
    margin-left: 0 !important; }
  .m-md-1 {
    margin: 0.25rem !important; }
  .mt-md-1,
  .my-md-1 {
    margin-top: 0.25rem !important; }
  .mr-md-1,
  .mx-md-1 {
    margin-right: 0.25rem !important; }
  .mb-md-1,
  .my-md-1 {
    margin-bottom: 0.25rem !important; }
  .ml-md-1,
  .mx-md-1 {
    margin-left: 0.25rem !important; }
  .m-md-2 {
    margin: 0.5rem !important; }
  .mt-md-2,
  .my-md-2 {
    margin-top: 0.5rem !important; }
  .mr-md-2,
  .mx-md-2 {
    margin-right: 0.5rem !important; }
  .mb-md-2,
  .my-md-2 {
    margin-bottom: 0.5rem !important; }
  .ml-md-2,
  .mx-md-2 {
    margin-left: 0.5rem !important; }
  .m-md-3 {
    margin: 1rem !important; }
  .mt-md-3,
  .my-md-3 {
    margin-top: 1rem !important; }
  .mr-md-3,
  .mx-md-3 {
    margin-right: 1rem !important; }
  .mb-md-3,
  .my-md-3 {
    margin-bottom: 1rem !important; }
  .ml-md-3,
  .mx-md-3 {
    margin-left: 1rem !important; }
  .m-md-4 {
    margin: 1.5rem !important; }
  .mt-md-4,
  .my-md-4 {
    margin-top: 1.5rem !important; }
  .mr-md-4,
  .mx-md-4 {
    margin-right: 1.5rem !important; }
  .mb-md-4,
  .my-md-4 {
    margin-bottom: 1.5rem !important; }
  .ml-md-4,
  .mx-md-4 {
    margin-left: 1.5rem !important; }
  .m-md-5 {
    margin: 3rem !important; }
  .mt-md-5,
  .my-md-5 {
    margin-top: 3rem !important; }
  .mr-md-5,
  .mx-md-5 {
    margin-right: 3rem !important; }
  .mb-md-5,
  .my-md-5 {
    margin-bottom: 3rem !important; }
  .ml-md-5,
  .mx-md-5 {
    margin-left: 3rem !important; }
  .p-md-0 {
    padding: 0 !important; }
  .pt-md-0,
  .py-md-0 {
    padding-top: 0 !important; }
  .pr-md-0,
  .px-md-0 {
    padding-right: 0 !important; }
  .pb-md-0,
  .py-md-0 {
    padding-bottom: 0 !important; }
  .pl-md-0,
  .px-md-0 {
    padding-left: 0 !important; }
  .p-md-1 {
    padding: 0.25rem !important; }
  .pt-md-1,
  .py-md-1 {
    padding-top: 0.25rem !important; }
  .pr-md-1,
  .px-md-1 {
    padding-right: 0.25rem !important; }
  .pb-md-1,
  .py-md-1 {
    padding-bottom: 0.25rem !important; }
  .pl-md-1,
  .px-md-1 {
    padding-left: 0.25rem !important; }
  .p-md-2 {
    padding: 0.5rem !important; }
  .pt-md-2,
  .py-md-2 {
    padding-top: 0.5rem !important; }
  .pr-md-2,
  .px-md-2 {
    padding-right: 0.5rem !important; }
  .pb-md-2,
  .py-md-2 {
    padding-bottom: 0.5rem !important; }
  .pl-md-2,
  .px-md-2 {
    padding-left: 0.5rem !important; }
  .p-md-3 {
    padding: 1rem !important; }
  .pt-md-3,
  .py-md-3 {
    padding-top: 1rem !important; }
  .pr-md-3,
  .px-md-3 {
    padding-right: 1rem !important; }
  .pb-md-3,
  .py-md-3 {
    padding-bottom: 1rem !important; }
  .pl-md-3,
  .px-md-3 {
    padding-left: 1rem !important; }
  .p-md-4 {
    padding: 1.5rem !important; }
  .pt-md-4,
  .py-md-4 {
    padding-top: 1.5rem !important; }
  .pr-md-4,
  .px-md-4 {
    padding-right: 1.5rem !important; }
  .pb-md-4,
  .py-md-4 {
    padding-bottom: 1.5rem !important; }
  .pl-md-4,
  .px-md-4 {
    padding-left: 1.5rem !important; }
  .p-md-5 {
    padding: 3rem !important; }
  .pt-md-5,
  .py-md-5 {
    padding-top: 3rem !important; }
  .pr-md-5,
  .px-md-5 {
    padding-right: 3rem !important; }
  .pb-md-5,
  .py-md-5 {
    padding-bottom: 3rem !important; }
  .pl-md-5,
  .px-md-5 {
    padding-left: 3rem !important; }
  .m-md-auto {
    margin: auto !important; }
  .mt-md-auto,
  .my-md-auto {
    margin-top: auto !important; }
  .mr-md-auto,
  .mx-md-auto {
    margin-right: auto !important; }
  .mb-md-auto,
  .my-md-auto {
    margin-bottom: auto !important; }
  .ml-md-auto,
  .mx-md-auto {
    margin-left: auto !important; } }
@media (min-width: 992px) {
  .m-lg-0 {
    margin: 0 !important; }
  .mt-lg-0,
  .my-lg-0 {
    margin-top: 0 !important; }
  .mr-lg-0,
  .mx-lg-0 {
    margin-right: 0 !important; }
  .mb-lg-0,
  .my-lg-0 {
    margin-bottom: 0 !important; }
  .ml-lg-0,
  .mx-lg-0 {
    margin-left: 0 !important; }
  .m-lg-1 {
    margin: 0.25rem !important; }
  .mt-lg-1,
  .my-lg-1 {
    margin-top: 0.25rem !important; }
  .mr-lg-1,
  .mx-lg-1 {
    margin-right: 0.25rem !important; }
  .mb-lg-1,
  .my-lg-1 {
    margin-bottom: 0.25rem !important; }
  .ml-lg-1,
  .mx-lg-1 {
    margin-left: 0.25rem !important; }
  .m-lg-2 {
    margin: 0.5rem !important; }
  .mt-lg-2,
  .my-lg-2 {
    margin-top: 0.5rem !important; }
  .mr-lg-2,
  .mx-lg-2 {
    margin-right: 0.5rem !important; }
  .mb-lg-2,
  .my-lg-2 {
    margin-bottom: 0.5rem !important; }
  .ml-lg-2,
  .mx-lg-2 {
    margin-left: 0.5rem !important; }
  .m-lg-3 {
    margin: 1rem !important; }
  .mt-lg-3,
  .my-lg-3 {
    margin-top: 1rem !important; }
  .mr-lg-3,
  .mx-lg-3 {
    margin-right: 1rem !important; }
  .mb-lg-3,
  .my-lg-3 {
    margin-bottom: 1rem !important; }
  .ml-lg-3,
  .mx-lg-3 {
    margin-left: 1rem !important; }
  .m-lg-4 {
    margin: 1.5rem !important; }
  .mt-lg-4,
  .my-lg-4 {
    margin-top: 1.5rem !important; }
  .mr-lg-4,
  .mx-lg-4 {
    margin-right: 1.5rem !important; }
  .mb-lg-4,
  .my-lg-4 {
    margin-bottom: 1.5rem !important; }
  .ml-lg-4,
  .mx-lg-4 {
    margin-left: 1.5rem !important; }
  .m-lg-5 {
    margin: 3rem !important; }
  .mt-lg-5,
  .my-lg-5 {
    margin-top: 3rem !important; }
  .mr-lg-5,
  .mx-lg-5 {
    margin-right: 3rem !important; }
  .mb-lg-5,
  .my-lg-5 {
    margin-bottom: 3rem !important; }
  .ml-lg-5,
  .mx-lg-5 {
    margin-left: 3rem !important; }
  .p-lg-0 {
    padding: 0 !important; }
  .pt-lg-0,
  .py-lg-0 {
    padding-top: 0 !important; }
  .pr-lg-0,
  .px-lg-0 {
    padding-right: 0 !important; }
  .pb-lg-0,
  .py-lg-0 {
    padding-bottom: 0 !important; }
  .pl-lg-0,
  .px-lg-0 {
    padding-left: 0 !important; }
  .p-lg-1 {
    padding: 0.25rem !important; }
  .pt-lg-1,
  .py-lg-1 {
    padding-top: 0.25rem !important; }
  .pr-lg-1,
  .px-lg-1 {
    padding-right: 0.25rem !important; }
  .pb-lg-1,
  .py-lg-1 {
    padding-bottom: 0.25rem !important; }
  .pl-lg-1,
  .px-lg-1 {
    padding-left: 0.25rem !important; }
  .p-lg-2 {
    padding: 0.5rem !important; }
  .pt-lg-2,
  .py-lg-2 {
    padding-top: 0.5rem !important; }
  .pr-lg-2,
  .px-lg-2 {
    padding-right: 0.5rem !important; }
  .pb-lg-2,
  .py-lg-2 {
    padding-bottom: 0.5rem !important; }
  .pl-lg-2,
  .px-lg-2 {
    padding-left: 0.5rem !important; }
  .p-lg-3 {
    padding: 1rem !important; }
  .pt-lg-3,
  .py-lg-3 {
    padding-top: 1rem !important; }
  .pr-lg-3,
  .px-lg-3 {
    padding-right: 1rem !important; }
  .pb-lg-3,
  .py-lg-3 {
    padding-bottom: 1rem !important; }
  .pl-lg-3,
  .px-lg-3 {
    padding-left: 1rem !important; }
  .p-lg-4 {
    padding: 1.5rem !important; }
  .pt-lg-4,
  .py-lg-4 {
    padding-top: 1.5rem !important; }
  .pr-lg-4,
  .px-lg-4 {
    padding-right: 1.5rem !important; }
  .pb-lg-4,
  .py-lg-4 {
    padding-bottom: 1.5rem !important; }
  .pl-lg-4,
  .px-lg-4 {
    padding-left: 1.5rem !important; }
  .p-lg-5 {
    padding: 3rem !important; }
  .pt-lg-5,
  .py-lg-5 {
    padding-top: 3rem !important; }
  .pr-lg-5,
  .px-lg-5 {
    padding-right: 3rem !important; }
  .pb-lg-5,
  .py-lg-5 {
    padding-bottom: 3rem !important; }
  .pl-lg-5,
  .px-lg-5 {
    padding-left: 3rem !important; }
  .m-lg-auto {
    margin: auto !important; }
  .mt-lg-auto,
  .my-lg-auto {
    margin-top: auto !important; }
  .mr-lg-auto,
  .mx-lg-auto {
    margin-right: auto !important; }
  .mb-lg-auto,
  .my-lg-auto {
    margin-bottom: auto !important; }
  .ml-lg-auto,
  .mx-lg-auto {
    margin-left: auto !important; } }
@media (min-width: 1200px) {
  .m-xl-0 {
    margin: 0 !important; }
  .mt-xl-0,
  .my-xl-0 {
    margin-top: 0 !important; }
  .mr-xl-0,
  .mx-xl-0 {
    margin-right: 0 !important; }
  .mb-xl-0,
  .my-xl-0 {
    margin-bottom: 0 !important; }
  .ml-xl-0,
  .mx-xl-0 {
    margin-left: 0 !important; }
  .m-xl-1 {
    margin: 0.25rem !important; }
  .mt-xl-1,
  .my-xl-1 {
    margin-top: 0.25rem !important; }
  .mr-xl-1,
  .mx-xl-1 {
    margin-right: 0.25rem !important; }
  .mb-xl-1,
  .my-xl-1 {
    margin-bottom: 0.25rem !important; }
  .ml-xl-1,
  .mx-xl-1 {
    margin-left: 0.25rem !important; }
  .m-xl-2 {
    margin: 0.5rem !important; }
  .mt-xl-2,
  .my-xl-2 {
    margin-top: 0.5rem !important; }
  .mr-xl-2,
  .mx-xl-2 {
    margin-right: 0.5rem !important; }
  .mb-xl-2,
  .my-xl-2 {
    margin-bottom: 0.5rem !important; }
  .ml-xl-2,
  .mx-xl-2 {
    margin-left: 0.5rem !important; }
  .m-xl-3 {
    margin: 1rem !important; }
  .mt-xl-3,
  .my-xl-3 {
    margin-top: 1rem !important; }
  .mr-xl-3,
  .mx-xl-3 {
    margin-right: 1rem !important; }
  .mb-xl-3,
  .my-xl-3 {
    margin-bottom: 1rem !important; }
  .ml-xl-3,
  .mx-xl-3 {
    margin-left: 1rem !important; }
  .m-xl-4 {
    margin: 1.5rem !important; }
  .mt-xl-4,
  .my-xl-4 {
    margin-top: 1.5rem !important; }
  .mr-xl-4,
  .mx-xl-4 {
    margin-right: 1.5rem !important; }
  .mb-xl-4,
  .my-xl-4 {
    margin-bottom: 1.5rem !important; }
  .ml-xl-4,
  .mx-xl-4 {
    margin-left: 1.5rem !important; }
  .m-xl-5 {
    margin: 3rem !important; }
  .mt-xl-5,
  .my-xl-5 {
    margin-top: 3rem !important; }
  .mr-xl-5,
  .mx-xl-5 {
    margin-right: 3rem !important; }
  .mb-xl-5,
  .my-xl-5 {
    margin-bottom: 3rem !important; }
  .ml-xl-5,
  .mx-xl-5 {
    margin-left: 3rem !important; }
  .p-xl-0 {
    padding: 0 !important; }
  .pt-xl-0,
  .py-xl-0 {
    padding-top: 0 !important; }
  .pr-xl-0,
  .px-xl-0 {
    padding-right: 0 !important; }
  .pb-xl-0,
  .py-xl-0 {
    padding-bottom: 0 !important; }
  .pl-xl-0,
  .px-xl-0 {
    padding-left: 0 !important; }
  .p-xl-1 {
    padding: 0.25rem !important; }
  .pt-xl-1,
  .py-xl-1 {
    padding-top: 0.25rem !important; }
  .pr-xl-1,
  .px-xl-1 {
    padding-right: 0.25rem !important; }
  .pb-xl-1,
  .py-xl-1 {
    padding-bottom: 0.25rem !important; }
  .pl-xl-1,
  .px-xl-1 {
    padding-left: 0.25rem !important; }
  .p-xl-2 {
    padding: 0.5rem !important; }
  .pt-xl-2,
  .py-xl-2 {
    padding-top: 0.5rem !important; }
  .pr-xl-2,
  .px-xl-2 {
    padding-right: 0.5rem !important; }
  .pb-xl-2,
  .py-xl-2 {
    padding-bottom: 0.5rem !important; }
  .pl-xl-2,
  .px-xl-2 {
    padding-left: 0.5rem !important; }
  .p-xl-3 {
    padding: 1rem !important; }
  .pt-xl-3,
  .py-xl-3 {
    padding-top: 1rem !important; }
  .pr-xl-3,
  .px-xl-3 {
    padding-right: 1rem !important; }
  .pb-xl-3,
  .py-xl-3 {
    padding-bottom: 1rem !important; }
  .pl-xl-3,
  .px-xl-3 {
    padding-left: 1rem !important; }
  .p-xl-4 {
    padding: 1.5rem !important; }
  .pt-xl-4,
  .py-xl-4 {
    padding-top: 1.5rem !important; }
  .pr-xl-4,
  .px-xl-4 {
    padding-right: 1.5rem !important; }
  .pb-xl-4,
  .py-xl-4 {
    padding-bottom: 1.5rem !important; }
  .pl-xl-4,
  .px-xl-4 {
    padding-left: 1.5rem !important; }
  .p-xl-5 {
    padding: 3rem !important; }
  .pt-xl-5,
  .py-xl-5 {
    padding-top: 3rem !important; }
  .pr-xl-5,
  .px-xl-5 {
    padding-right: 3rem !important; }
  .pb-xl-5,
  .py-xl-5 {
    padding-bottom: 3rem !important; }
  .pl-xl-5,
  .px-xl-5 {
    padding-left: 3rem !important; }
  .m-xl-auto {
    margin: auto !important; }
  .mt-xl-auto,
  .my-xl-auto {
    margin-top: auto !important; }
  .mr-xl-auto,
  .mx-xl-auto {
    margin-right: auto !important; }
  .mb-xl-auto,
  .my-xl-auto {
    margin-bottom: auto !important; }
  .ml-xl-auto,
  .mx-xl-auto {
    margin-left: auto !important; } }
.text-justify {
  text-align: justify !important; }
.text-nowrap {
  white-space: nowrap !important; }
.text-truncate {
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap; }
.text-left {
  text-align: left !important; }
.text-right {
  text-align: right !important; }
.text-center {
  text-align: center !important; }
@media (min-width: 576px) {
  .text-sm-left {
    text-align: left !important; }
  .text-sm-right {
    text-align: right !important; }
  .text-sm-center {
    text-align: center !important; } }
@media (min-width: 768px) {
  .text-md-left {
    text-align: left !important; }
  .text-md-right {
    text-align: right !important; }
  .text-md-center {
    text-align: center !important; } }
@media (min-width: 992px) {
  .text-lg-left {
    text-align: left !important; }
  .text-lg-right {
    text-align: right !important; }
  .text-lg-center {
    text-align: center !important; } }
@media (min-width: 1200px) {
  .text-xl-left {
    text-align: left !important; }
  .text-xl-right {
    text-align: right !important; }
  .text-xl-center {
    text-align: center !important; } }
.text-lowercase {
  text-transform: lowercase !important; }
.text-uppercase {
  text-transform: uppercase !important; }
.text-capitalize {
  text-transform: capitalize !important; }
.font-weight-light {
  font-weight: 300 !important; }
.font-weight-normal {
  font-weight: 400 !important; }
.font-weight-bold {
  font-weight: 700 !important; }
.font-italic {
  font-style: italic !important; }
.text-white {
  color: #fff !important; }
.text-primary {
  color: #512479 !important; }
a.text-primary:hover, a.text-primary:focus {
  color: #49075e !important; }
.text-secondary {
  color: #333 !important; }
a.text-secondary:hover, a.text-secondary:focus {
  color: #1a1a1a !important; }
.text-success {
  color: #28a745 !important; }
a.text-success:hover, a.text-success:focus {
  color: #1e7e34 !important; }
.text-info {
  color: #17a2b8 !important; }
a.text-info:hover, a.text-info:focus {
  color: #117a8b !important; }
.text-warning {
  color: #ffc107 !important; }
a.text-warning:hover, a.text-warning:focus {
  color: #d39e00 !important; }
.text-danger {
  color: #dc3545 !important; }
a.text-danger:hover, a.text-danger:focus {
  color: #bd2130 !important; }
.text-light {
  color: #f8f9fa !important; }
a.text-light:hover, a.text-light:focus {
  color: #dae0e5 !important; }
.text-dark {
  color: #343a40 !important; }
a.text-dark:hover, a.text-dark:focus {
  color: #1d2124 !important; }
.text-muted {
  color: #6c757d !important; }
.text-hide {
  font: 0/0 a;
  color: transparent;
  text-shadow: none;
  background-color: transparent;
  border: 0; }
.visible {
  visibility: visible !important; }
.invisible {
  visibility: hidden !important; }
@media print {
  *,
  *::before,
  *::after {
    text-shadow: none !important;
    box-shadow: none !important; }
  a:not(.btn) {
    text-decoration: underline; }
  abbr[title]::after {
    content: " (" attr(title) ")"; }
  pre {
    white-space: pre-wrap !important; }
  pre,
  blockquote {
    border: 1px solid #999;
    page-break-inside: avoid; }
  thead {
    display: table-header-group; }
  tr,
  img {
    page-break-inside: avoid; }
  p,
  h2,
  h3 {
    orphans: 3;
    widows: 3; }
  h2,
  h3 {
    page-break-after: avoid; }
  @page {
    size: a3; }
  body {
    min-width: 992px !important; }
  .container {
    min-width: 992px !important; }
  .navbar {
    display: none; }
  .badge {
    border: 1px solid #000; }
  .table {
    border-collapse: collapse !important; }
    .table td,
    .table th {
      background-color: #fff !important; }
  .table-bordered th,
  .table-bordered td {
    border: 1px solid #ddd !important; } }
/********************************************************************
    ZYPOP - HTTPS://ZYPOPWEBTEMPLATES.COM
    FREE WEB TEMPLATES
********************************************************************/
/**
This file contains the core CSS for this template, built on the Bootstrap framework
Set the variable values in _variables.scss
**/
/**
Page structure and elements
**/
.bg-primary {
  background-color: #512479 !important; }
body {
  background: #512479; }
main {
  padding-top: 1em; }
.text-primary {
  color: #512479 !important; }
.container-fluid {
  background-color: #fff; }
#sidebar {
  width: auto;
  padding: 0 15px;
  position: auto;
  overflow: scroll;
  height: 100%; }
.right-sidebar #sidebar {
  right: 0; }
#content {
  margin-left: 0; }
a, .page-link {
  color: #512479, 10%;
  font-weight: bold; }
  a:hover, a:focus, a:active, .page-link:hover, .page-link:focus, .page-link:active {
    color: #49075e; }
.page-item.active .page-link {
  border-color: #512479;
  background-color: #512479; }
fieldset {
  margin-bottom: 1rem;
  display: block;
  border-top: 1px solid #ccc; }
fieldset legend {
  width: auto;
  padding-right: 0.5rem;
  font-size: 1.1rem;
  font-weight: bold; }
table th {
  background-color: #512479;
  color: #fff;
  border-color: #523047 !important;
  border-bottom: 1px solid #523047; }
blockquote {
  display: block;
  border-left: 5px solid #ccc;
  padding: 0.5rem;
  color: #666;
  margin-bottom: 1rem; }
.btn {
  font-weight: bold; }
.btn-secondary, .badge-secondary,
.btn-primary, .badge-primary {
  color: #fff; }
.btn-outline-primary:hover, .btn-outline-primary:focus, .btn-outline-primary:active {
  color: #fff; }
h1, h2, h3, h4, h5, h6, .h1, .h2, .h3, .h4, .h5, .h6 {
  color: #333;
  font-weight: bold; }
/**
Jumbotron / Slider
**/
.jumbotron-narrow {
  padding: 2rem; }
.jumbotron-wrap .container {
  padding-bottom: 1rem; }
#mainCarousel, .static-slider {
  border-bottom: 5px solid #ebebeb;
  border-radius: 0px; }
.jumbotron-wrap .jumbotron {
  background: #f8f8f8;
  margin-bottom: 0; }
.jumbotron-wrap h1, .jumbotron-wrap .h1 {
  color: #5f5f5f; }
/**
Header
**/
header {
  margin-top: 3.5rem;
  padding: 0 30px 15px; }
header h1 {
  color: #2d1927;
  font-size: 2.5rem;
  text-align: left;
  line-height: 2.5rem;
  letter-spacing: -0.1rem;
  margin-bottom: 0; }
header h1 span {
  color: #fff; }
.sidebar-social-icons {
  margin-top: 20px; }
.sidebar-social-icons a {
  color: white; }
  .sidebar-social-icons a:hover, .sidebar-social-icons a:focus, .sidebar-social-icons a:active {
    color: #333; }
/**
Navbar
**/
.navbar-toggler {
  margin: 0.5rem; }
.navbar, #mainNavbar {
  padding: 0;
  width: 100%; }
.navbar ul {
  padding: 0;
  list-style: none; }
.navbar ul.sub-navbar {
  padding-left: 20px;
  background-color: #2d1927; }
.navbar .active ul.sub-navbar {
  background-color: #fff; }
.mobile-header-controls {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  justify-content: space-between; }
#mainNavbar .nav-link {
  padding: 15px 30px;
  color: #fff;
  background-color: #49075e; }
  #mainNavbar .nav-link:hover, #mainNavbar .nav-link:focus, #mainNavbar .nav-link:active {
    text-decoration: underline; }
#mainNavbar ul.sub-navbar .nav-link {
  background-color: #2d1927; }
#mainNavbar .active .nav-link {
  color: #512479;
  background-color: #fff !important; }
.navbar-dark .navbar-brand {
  color: #fff; }
  .navbar-dark .navbar-brand:hover, .navbar-dark .navbar-brand:focus, .navbar-dark .navbar-brand:active {
    color: #fff; }
.navbar-dark .navbar-brand span {
  color: #2d1927; }
/**
Footer
**/
.footer-container {
  padding-top: 1rem;
  padding-bottom: 1rem; }
footer {
  color: #acacac;
  font-size: 0.9rem; }
.footer-lists {
  background-color: #f8f8f8;
  padding: 30px;
  margin-bottom: 1rem;
  border-top: 5px solid #ebebeb; }
.footer-bottom {
  color: #c5c5c5; }
.footer-bottom a {
  color: #acacac; }
  .footer-bottom a:hover, .footer-bottom a:focus, .footer-bottom a:active {
    color: #9f9f9f;
    border-bottom-color: #9f9f9f; }
.footer-lists h4 {
  color: #929292; }
.social-icons a {
  margin-right: 15px;
  border-bottom: none; }
.footer-lists ul {
  list-style: none;
  margin: 0;
  padding: 0; }
.footer-lists ul li {
  padding: 0.2rem 0; }
footer p {
  margin: 0;
  padding-bottom: 1rem; }
footer p:last-child {
  padding-bottom: 0; }
footer a {
  color: #b8b8b8;
  border-bottom: 1px solid #b8b8b8;
  font-weight: normal; }
  footer a:hover, footer a:focus, footer a:active {
    text-decoration: none;
    border-bottom-color: #5f5f5f;
    color: #5f5f5f; }
/**
Articles
**/
article {
  margin-bottom: 2rem;
  border-bottom: 1px solid #dee2e6;
  padding-bottom: 2rem; }
article h2.article-title {
  font-size: 2.5rem;
  margin-bottom: 0;
  color: #333;
  letter-spacing: -1px; }
article p.article-meta {
  color: #ccc;
  font-size: 0.8rem; }
/**
Sidebar
**/
.sidebar-box {
  margin-bottom: 2rem; }
.sidebar-box-bg {
  padding: 1rem;
  background-color: #f8f8f8;
  border-radius: 0px; }
/***
Better nesting of list groups
***/
.list-group-item {
  border: none;
  border-bottom: 3px solid #fff; }
.sidebar-box-bg a, .list-group-item {
  color: #5f5f5f; }
  .sidebar-box-bg a:hover, .sidebar-box-bg a:focus, .sidebar-box-bg a:active, .list-group-item:hover, .list-group-item:focus, .list-group-item:active {
    color: #464646;
    text-decoration: underline; }
.list-group .list-group .list-group-item {
  padding-left: 2.5rem; }
.list-group .list-group .list-group .list-group-item {
  padding-left: 3.75rem; }
.list-group .list-group .list-group .list-group .list-group-item {
  padding-left: 5rem; }
.list-group > .list-group .list-group-item:first-child,
.list-group > .list-group .list-group-item:last-child {
  border-radius: 0; }
.list-group > .list-group .list-group-item:last-child {
  border-bottom: 2px solid #fff; }
.list-group-root {
  background-color: #f8f8f8;
  padding: 0rem;
  border-radius: 0px; }
.list-group-item {
  background-color: #f8f8f8; }
.list-group-item.active {
  background-color: #512479;
  border-color: #fff;
  color: #fff;
  border-radius: 0px !important; }
  .list-group-item.active:hover, .list-group-item.active:focus, .list-group-item.active:active {
    color: #fff;
    text-decoration: underline; }
/**
Responsive typography
https://getbootstrap.com/docs/4.0/content/typography/#responsive-typography
**/
html {
  font-size: 16px; }
.navbar-container {
  padding: 0; }
html {
  font-size: 12px; }
@media (min-width: 768px) {
  html {
    font-size: 14px; }
  #sidebar {
    width: 320px;
    position: fixed;
    padding: 0; }
  #content {
    padding-left: 320px; }
  .right-sidebar #content {
    padding-left: 0 !important;
    padding-right: 320px; }
  #content #content-wrapper {
    background-color: #fff;
    padding: 15px; } }
@media (min-width: 992px) {
  html {
    font-size: 16px; } }''')


def create_header(bin_list, directory, active, long_read_qc_html, short_read_qc_html=None):
    main_header = '''<!doctype html>
<html lang="en">
    <head>
        <title>slamM</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <!-- Main CSS -->
        <link rel="stylesheet" href="css/style.css">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/v/bs4/jq-3.3.1/jszip-2.5.0/dt-1.10.18/b-1.5.6/b-html5-1.5.6/cr-1.5.0/datatables.min.css"/>
        <script src="https://code.jquery.com/jquery-3.3.1.min.js" integrity="sha256-FgpCb/KJQlLNfOu91ta32o/NMZxltwRo8QtmkMRdAu8=" crossorigin="anonymous"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
        <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/pdfmake.min.js"></script>
        <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/vfs_fonts.js"></script>
        <script type="text/javascript" src="https://cdn.datatables.net/v/bs4/jq-3.3.1/jszip-2.5.0/dt-1.10.18/b-1.5.6/b-html5-1.5.6/cr-1.5.0/datatables.min.js"></script>
    </head>
    <body>


        <!-- Main navigation -->
        <div id="sidebar">

            <div class="navbar-expand-md navbar-dark"> 

                <header class="d-none d-md-block">
                    <h1><span>slam</span>M</h1>
                </header>


                <!-- Mobile menu toggle and header -->
                <div class="mobile-header-controls">
                    <a class="navbar-brand d-md-none d-lg-none d-xl-none" href="#"><span>slam</span>M</a>
                    <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#SidebarContent" aria-controls="SidebarContent" aria-expanded="false" aria-label="Toggle navigation">
                        <span class="navbar-toggler-icon"></span>
                    </button>
                </div>

                <div id="SidebarContent" class="collapse flex-column navbar-collapse">



                    <!-- Main navigation items -->
                    <nav class="navbar navbar-dark">
                        <div id="mainNavbar">
                            <ul class="flex-column mr-auto">
'''
    if active == 'index.html':
        main_header += '                                <li class="nav-item active">\n'
    else:
        main_header += '                                <li class="nav-item">\n'

    main_header += '''                                    <a class="nav-link" href="''' + directory + '''index.html">Home <span class="sr-only">(current)</span></a>
                                </li>
                                <li class="nav-item">
                                        <a class="nav-link" href="''' + directory + long_read_qc_html + '''">Long read stats</a>
                                </li>'''

    if not short_read_qc_html is None:
        main_header += '''                                <li class="nav-item">
                                        <a class="nav-link" href="''' + directory + short_read_qc_html + '''">Short read stats</a>
                                </li>
'''
    main_header += '''                                <li class="nav-item dropdown">
                                            <a class="nav-link dropdown-toggle" href="#MenuDropdown" data-toggle="collapse" aria-controls="MenuDropdown" aria-expanded="false">Details of bins</a>
                                            <ul id="MenuDropdown" class="sub-navbar collapse flex-column">'''
    for i in bin_list:
        if active == i:
            main_header += '                                                <li class="nav-item active"><a class="nav-link" href="%sbin/%s.html">Bin %s</a></li>\n' % (
            directory, i, i)
        else:
            main_header += '                                                <li class="nav-item"><a class="nav-link" href="%sbin/%s.html">Bin %s</a></li>\n' % (
            directory, i, i)
    main_header += '''                                            </ul>
                                </li>
                                <li class="nav-item">
                                        <a class="nav-link" href="gtdbtk.html">GTDBtk</a>
                                </li>
                                <li class="nav-item">
                                        <a class="nav-link" href="https://github.com/mjsull/SDMass/issues">Help</a>
                                </li>
                            </ul>
                        </div>   
                    </nav>

                </div>
            </div> 
        </div>        


        <div id="content">
            <div id="content-wrapper">
'''

    return (main_header)


def add_title(header, subheader=''):
    return ('''
            <!-- Jumbtron / Slider -->
                <div class="jumbotron-wrap">
                    <div class="container-fluid">
                        <div class="jumbotron static-slider">
                            <h1 class="text-center">''' + header + \
            '''</h1>
                                   <p class="lead text-center">''' + subheader + '''</p>
                        </div>
                    </div>
                </div>
''')


def add_main(header, text, contigs=None):
    if contigs is None:
        return ('''                <main class="container-fluid">
                    <div class="row">
                        <!-- Main content -->
                        <div class="col-md-12">
                            <article>
                                <h2 class="article-title">''' + header + '''</h2>
                                <p> ''' + text + ''' </p>
''')
    else:
        html_string = '''                <main class="container-fluid">
                    <div class="row">
                        <!-- Sidebar -->
                        <aside class="col-md-2">
                            <div class="sidebar-box">
                                <h4>Contigs</h4>
                                <div class="list-group list-group-root">
                                    <a class="list-group-item active" href="index.html">Overview</a>
'''
        count = 0
        for i in contigs:
            the_path = "www/contigs/" + i + ".html"
            if os.path.exists(the_path):
                html_string += '                                    <a class="list-group-item" href="' + the_path + '">' + i + '</a>\n'
            else:
                count += 1
        html_string += '''                                </div>
                            </div>
                            <div class="sidebar-box sidebar-box-bg">
                                <h4>n.b.</h4>
                                <p> ''' + str(count) + ''' contigs under 100Kbp not shown.</p>
                            </div>
                        </aside>
                        <!-- Main content -->
                        <div class="col-md-10">
                            <article>
                                <h2 class="article-title">''' + header + '''</h2>
                                <p> ''' + text + ''' </p>
'''
        return (html_string)


def end_main():
    return ('''
                        </div>
                    </div> 
                </main>
''')


def add_footer():
    return ('''                <!-- Footer -->
                <div class="container-fluid footer-container">
                    <footer class="footer">
                        <div class="footer-bottom">
                                <p class="text-center">Created by <a href="https://mjsull.github.io"</a>mjsull.github.io</p>
                                <p class="text-center"><a href="#">Back to top</a></p>
                        </div>
                    </footer>
                </div> 
            </div>
        </div>
        <!-- Bootcamp JavaScript -->
        <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    </body>
</html>''')


def create_table(headers, list_of_vals, hide_by_default=set(), text_table=None):
    if not text_table is None:
        with open(text_table, 'w') as o:
            o.write('\t'.join(headers) + '\n')
            for i in list_of_vals:
                o.write('\t'.join(map(str, i)) + '\n')
    out_string = '''<script>
    $(document).ready(function() {
    var table = $('#thetable').DataTable( {
        dom: 'frtiplB;',
        "scrollX": true,
        buttons: [
            'copy', 'csv', 'excel', 'pdf'
        ]
    });

    $('a.toggle-vis').on( 'click', function (e) {
        e.preventDefault();
        $(this).toggleClass('greened');
        // Get the column API object
        var column = table.column( $(this).attr('data-column') );

        // Toggle the visibility
        column.visible( ! column.visible() );
    } );
'''
    for num, head in enumerate(headers):
        if head in hide_by_default:
            out_string += '        table.column(%d).visible( false );\n' % (num)
    out_string += '''
} );
    </script>
    <div>
        Toggle column:
'''
    for num, head in enumerate(headers):
        if head in hide_by_default:
            out_string += '<a class="toggle-vis greened" data-column="%d">%s</a> - ' % (num, head)
        else:
            out_string += '<a class="toggle-vis" data-column="%d">%s</a> - ' % (num, head)
    out_string = out_string[:-3]
    out_string += '\n        </div>'
    out_string += '''                                <table class="table" id="thetable" style="width:100%">
                                    <thead>
                                        <tr>
'''
    for i in headers:
        out_string += '                                            <th>' + i + '</th>\n'
    out_string += '''                                        </tr>
                                    </thead>
                                    <tbody>
'''
    for i in list_of_vals:
        out_string += '                                        <tr>\n'
        for j in i:
            out_string += '                                            <td>' + str(j) + '</td>\n'
        out_string += '                                        </tr>\n'
    out_string += '''                                    </tbody>
                                </table>
'''
    return (out_string)


def get_cov_stats_long(bamfile, contig, bin_size=3000, bin_step=500, buffer=50):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    ref_length = samfile.get_reference_length(contig)
    bin_num = ref_length // bin_step + 1
    coverage_forward = numpy.zeros(bin_num)
    coverage_reverse = numpy.zeros(bin_num)
    trimmed_starts = numpy.zeros(bin_num, dtype=int)
    trimmed_ends = numpy.zeros(bin_num, dtype=int)
    through = numpy.zeros(bin_num, dtype=int)
    starts_in = numpy.zeros(bin_num, dtype=int)
    ends_in = numpy.zeros(bin_num, dtype=int)
    x = [y for y in range(0, ref_length + 1, bin_num)]
    for read in samfile.fetch(contig):
        rstart = read.reference_start
        rend = read.reference_end
        qstart = read.query_alignment_start
        qend = read.query_alignment_end
        qlength = read.infer_query_length()
        if qstart > buffer:
            trimmed_start = True
        else:
            trimmed_start = False
        if qend < qlength - buffer:
            trimmed_end = True
        else:
            trimmed_end = False
        if read.is_secondary:
            continue
        for i in range(max([0, (rstart - bin_size) // bin_step]), rend // bin_step + 1):
            if i * bin_step <= rstart < rend < i * bin_step + bin_size:
                pass
            elif i * bin_step <= rstart < i * bin_step + bin_size:
                if trimmed_start:
                    trimmed_starts[i] += 1
                else:
                    starts_in[i] += 1
            elif i * bin_step <= rend < i * bin_step + bin_size:
                if trimmed_end:
                    trimmed_ends[i] += 1
                else:
                    ends_in[i] += 1
            elif rstart <= i * bin_step and rend >= i * bin_step + bin_size:
                through[i] += 1
            bases = (min([rend, i * bin_step + bin_size]) - max([rstart, i * bin_step])) / bin_size
            if read.is_reverse:
                coverage_reverse[i] += bases
            else:
                coverage_forward[i] += bases
    return (coverage_forward, coverage_reverse, trimmed_starts, trimmed_ends, starts_in, ends_in, x)


def get_gtdbtk(gtdbtk_folder, in_dict=None):
    connect_dict = {}
    out_dict = {}
    bins = 0
    # get bacterial information
    if os.path.exists(os.path.join(gtdbtk_folder, 'gtdbtk.bac120.summary.tsv')):
        with open(os.path.join(gtdbtk_folder, 'gtdbtk.bac120.summary.tsv')) as f:
            f.readline()
            for line in f:
                bins += 1
                the_bin, phylo, nearest, ani_radius, ani_tax, ani = line.split('\t')[:6]
                # the_bin = the_bin.split('.')[-1]
                out_dict[the_bin] = (nearest, ani, phylo)
                if not in_dict is None:
                    cov = in_dict[the_bin]
                else:
                    cov = 10
                lastname = None
                # add bin to dictionary of connections and add weight of bin to all upstream connections
                for i in phylo.split(';'):
                    if len(i) == 3:
                        the_name = i + the_bin
                    else:
                        the_name = i
                    if not lastname is None:
                        if (lastname, the_name) in connect_dict:
                            connect_dict[(lastname, the_name)] += cov
                        else:
                            connect_dict[(lastname, the_name)] = cov
                    lastname = the_name
    # do the same for archaea
    if os.path.exists(os.path.join(gtdbtk_folder, 'gtdbtk.ar122.summary.tsv')):
        with open(os.path.join(gtdbtk_folder, 'gtdbtk.ar122.summary.tsv')) as f:
            f.readline()
            for line in f:
                bins += 1
                the_bin, phylo, nearest, ani_radius, ani_tax, ani = line.split('\t')[:6]
                #  the_bin = the_bin.split('.')[-1]
                out_dict[the_bin] = (nearest, ani, phylo)
                if not in_dict is None:
                    cov = in_dict[the_bin]
                else:
                    cov = 10
                lastname = None
                for i in phylo.split(';'):
                    if len(i) == 3:
                        the_name = i + the_bin
                    else:
                        the_name = i
                    if not lastname is None:
                        if (lastname, the_name) in connect_dict:
                            connect_dict[(lastname, the_name)] += cov
                        else:
                            connect_dict[(lastname, the_name)] = cov
                    lastname = the_name
    # if not provided with weights just return a dictionary of information
    if in_dict is None:
        return out_dict
    connect_list = []
    # turn the dictionary into a list of connections
    for i in connect_dict:
        connect_list.append([i[0], i[1], connect_dict[i]])
    # quit gracefully if there's no output from gtdbtk
    if connect_list == []:
        with open('www/gtdbtk.html', 'w') as o:
            o.write("No gtdbtk output")
        return
    depth_first = []
    # sort list by weight
    connect_list.sort(key=lambda x: x[2], reverse=True)
    todo = ['d__Archaea', 'd__Bacteria']
    # order list with a depth first search
    while todo != []:
        getit = todo.pop()
        to_add = []
        for i in connect_list:
            if i[0] == getit:
                depth_first.append(i)
                to_add.append(i[1])
        for i in to_add[::-1]:
            todo.append(i)
    connect_list = depth_first
    # create the webpage with sankey diagram, iterations set to 0 to preserve order provided to the webpage
    with open('www/gtdbtk.html', 'w') as o:
        o.write('''<html>
<body>
 <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
<div id="sankey_multiple" style="width: 800px; height: 3000px;"></div>
<script type="text/javascript">
  google.charts.load("current", {packages:["sankey"]});
  google.charts.setOnLoadCallback(drawChart);
   function drawChart() {
    var data = new google.visualization.DataTable();
    data.addColumn('string', 'From');
    data.addColumn('string', 'To');
    data.addColumn('number', 'Weight');
    data.addRows([\n''')
        for i in connect_list[:-1]:
            o.write(str(i) + ',\n')
        o.write(str(connect_list[-1]))
        o.write('''
    ]);
    // Set chart options
    var options = {
      width: 1400,
      height: ''' + str(bins * 30) +
                ''', 
                      sankey: {
                          iterations: 0,
                          },
                      textStyle: {
                            fontSize: 6,
                        },
                    };
                    // Instantiate and draw our chart, passing in some options.
                    var chart = new google.visualization.Sankey(document.getElementById('sankey_multiple'));
                    chart.draw(data, options);
                   }
                </script>
                </body>
                </html>''')


def get_cov_stats_short(bamfile, contig, bin_size=3000, bin_step=500):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    ref_length = samfile.get_reference_length(contig)
    bin_num = ref_length // bin_step + 1
    coverage_forward = numpy.zeros(bin_num)
    coverage_reverse = numpy.zeros(bin_num)
    for read in samfile.fetch(contig):
        rstart = read.reference_start
        rend = read.reference_end
        if read.is_secondary or rstart is None or rend is None:
            continue
        for i in range(max([0, (rstart - bin_size) // bin_step]), rstart // bin_step + 1):
            bases = (min([rend, i * bin_step + bin_size]) - max([rstart, i * bin_step])) / bin_size
            if read.is_reverse:
                coverage_reverse[i] += bases
            else:
                coverage_forward[i] += bases
    return coverage_forward, coverage_reverse


def get_gene_sizes(gff_file):
    size_dict = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            contig, prod, feat, start, stop = line.split()[:5]
            if not contig in size_dict:
                size_dict[contig] = []
            size_dict[contig].append(int(stop) - int(start))
    return size_dict


# def create_contig_page(bin, ctg, cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in, ends_ind,
#                        cov_forward_ill, cov_reverse_ill, x, outpage):
#      html_out = open(out_directory + '/qc_website/ctgs/' + str(html_name) + '_graphs.html', 'w')
#      html_out.write(header)
#      html_out.write('  <script type="text/javascript">\n'
#                     '  window.onload = function () {\n'
#                     '  var chart1 = new CanvasJS.Chart("chartContainer1",\n'
#                    '    {\n'
#                    '      zoomEnabled: true,\n'
#                    '      title:{\n'
#                    '      text: "Total number of each read type in bin - bin size: ' + str(
#         bin_size) + ', bin step: ' + str(bin_step) + ')",\n'
#                                                      '      fontSize: 24\n'
#                                                      '      },\n'
#                                                      '      axisX: {\n'
#                                                      '      title: "Position in genome",\n'
#                                                      '      titleFontSize: 16,\n')
# #         if clipped_flag != [] or unclipped_flag != []:
# #             html_out.write('      stripLines:[\n')
# #             for values in clipped_flag[:-1]:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(values[0]) + ',\n'
# #                                                                          '        endValue: ' + str(
# #                     values[1]) + ',\n'
# #                                  '        color: "#16CBEF"\n'
# #                                  '      },\n')
# #             if clipped_flag != []:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(clipped_flag[-1][0]) + ',\n'
# #                                                                                    '        endValue: ' + str(
# #                     clipped_flag[-1][1]) + ',\n'
# #                                            '        color: "#16CBEF"\n'
# #                                            '      }')
# #                 if unclipped_flag != []:
# #                     html_out.write(',')
# #                 else:
# #                     html_out.write('\n')
# #             for values in unclipped_flag[:-1]:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(values[0]) + ',\n'
# #                                                                          '        endValue: ' + str(
# #                     values[1]) + ',\n'
# #                                  '        color: "#53A01F"\n'
# #                                  '      },\n')
# #             if unclipped_flag != []:
# #                 html_out.write('      {\n'
# #                                '        startValue: ' + str(unclipped_flag[-1][0]) + ',\n'
# #                                                                                      '        endValue: ' + str(
# #                     unclipped_flag[-1][1]) + ',\n'
# #                                              '        color: "#53A01F"\n'
# #                                              '      }\n')
# #             html_out.write('       ]\n  ')
# #
#     html_out.write('      },\n'
#                    '      axisY: {\n'
#                    '      title: "Number of reads",\n'
#                    '      titleFontSize: 16\n'
#                    '      },\n'
#                    '       data: [')
#     for the_data, leg_lab, colour in zip(
#             [forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end,
#              forward_start_clipped,
#              reverse_start_clipped, forward_end_clipped, reverse_end_clipped], ['Read spans bin (Forward)',
#                                                                                 'Read spans bin (Reverse)',
#                                                                                 'Read starts in bin (F)',
#                                                                                 'Read starts in bin (R)',
#                                                                                 'Read terminates in bin (F)',
#                                                                                 'Read terminates in bin (R)',
#                                                                                 'Read start clipped in bin (F)',
#                                                                                 'Read start clipped in bin (R)',
#                                                                                 'Read end clipped in bin (F)',
#                                                                                 'Read end clipped in bin (R)'],
#             ["#BE4200", "#EE16B8", "#38081F", "#4282F6", "#45FAA5", "#763181",
#              "#234D65", "#F1F272", "#D15AF4", "#CBC22F"]):
#         html_out.write('\n{\n'
#                        '         type: "line",\n'
#                        '         showInLegend: true,\n'
#                        '         legendText: "' + leg_lab + '",\n'
#                                                             '         color: "' + colour + '",\n'
#                                                                                            '         markerType: "none",\n'
#                                                                                            '         dataPoints: [\n')
#         if refnum == 0:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'w')
#             bgt_out = open(
#                 options.output_folder + '/bigwig/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.bwt',
#                 'w')
#             bgt_out.write(
#                 'track type=bigwig bigDataUrl=https://vanbah01.u.hpc.mssm.edu/igb/' + options.assembly_name + '/'
#                 + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower()
#                 + '.bw name=' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() +
#                 'color=0,0,200 altColor=0,200,0 autoScale=on alwaysZero=on graphType=bar yLineMark=10 yLineOnOff=on\n')
#             bgt_out.close()
#         else:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'a')
#         bg_out.write(
#             'fixedStep chrom=' + reference + ' start=1 step=' + str(bin_step) + ' span=' + str(bin_step) + '\n')
#         for value in range(0, len(x_axis) - 1):
#             html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
#             if x_axis[value] >= 0 and x_axis[value] < chrom_size[reference] - bin_step:
#                 bg_out.write(str(the_data[value]) + '\n')
#         bg_out.close()
#         html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
#         html_out.write('        ]\n      }')
#         if not the_data is reverse_end_clipped:
#             html_out.write(',\n')
#     html_out.write('      ],\n'
#                    '      rangeChanged: syncHandler\n'
#                    '    });\n'
#                    '    chart1.render();\n'
#                    '  var chart2 = new CanvasJS.Chart("chartContainer2",\n'
#                    '    {\n'
#                    '      zoomEnabled: true,\n'
#                    '      title:{\n'
#                    '      text: "Proportion of each read type in bin - bin size: ' + str(
#         bin_size) + ', bin step: ' + str(bin_step) + ')",\n'
#                                                      '      fontSize: 24\n'
#                                                      '      },\n'
#                                                      '      axisX: {\n'
#                                                      '      title: "Position in genome",\n'
#                                                      '      titleFontSize: 16\n'
#                                                      '},\n'
#                                                      '      axisY: {\n'
#                                                      '      title: "Number of reads",\n'
#                                                      '      titleFontSize: 16\n'
#                                                      '},\n'
#                                                      '       data: [')
#     for the_data, leg_lab in zip(
#             [forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end,
#              forward_start_clipped,
#              reverse_start_clipped, forward_end_clipped, reverse_end_clipped], ['Read spans bin (Forward)',
#                                                                                 'Read spans bin (Reverse)',
#                                                                                 'Read starts in bin (F)',
#                                                                                 'Read starts in bin (R)',
#                                                                                 'Read terminates in bin (F)',
#                                                                                 'Read terminates in bin (R)',
#                                                                                 'Read start clipped in bin (F)',
#                                                                                 'Read start* clipped in bin (R)',
#                                                                                 'Read end clipped in bin (F)',
#                                                                                 'Read end clipped in bin (R)']):
#         html_out.write('         {\n'
#                        '         type: "stackedArea100",\n'
#                        '         showInLegend: true,\n'
#                        '         legendText: "' + leg_lab + '",\n'
#                                                             '         markerType: "none",\n'
#                                                             '         legendMarkerType: "square",\n'
#                                                             '        dataPoints: [\n')
#         for value in range(0, len(x_axis) - 1):
#             html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
#
#         html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
#         html_out.write('        ]\n      }')
#         if not the_data is reverse_end_clipped:
#             html_out.write(',')
#     html_out.write('      ],\n'
#                    '      rangeChanged: syncHandler\n'
#                    '   });\n'
#                    '    chart2.render();\n'
#                    '  var chart3 = new CanvasJS.Chart("chartContainer3",\n'
#                    '    {\n'
#                    '      zoomEnabled: true,\n'
#                    '      title:{\n'
#                    '      text: "Number of large indels in each bin - bin size: ' + str(
#         bin_size) + ', bin step: ' + str(bin_step) + '",\n'
#                                                      '      fontSize: 24\n'
#                                                      '      },\n'
#                                                      '      axisX: {\n'
#                                                      '      title: "Position in genome",\n'
#                                                      '      titleFontSize: 16,\n')
#     if indel_flag != []:
#         html_out.write('      stripLines:[\n')
#         for values in indel_flag[:-1]:
#             html_out.write('      {\n'
#                            '        startValue: ' + str(values[0]) + ',\n'
#                                                                      '        endValue: ' + str(
#                 values[1]) + ',\n'
#                              '        color: "#16CBEF"\n'
#                              '      },\n')
#         html_out.write('      {\n'
#                        '        startValue: ' + str(indel_flag[-1][0]) + ',\n'
#                                                                          '        endValue: ' + str(
#             indel_flag[-1][1]) + ',\n'
#                                  '        color: "#16CBEF"\n'
#                                  '      }\n'
#                                  '      ]\n')
#     html_out.write('\n    },\n'
#                    '      axisY: {\n'
#                    '      title: "Number of reads",\n'
#                    '      titleFontSize: 16\n'
#                    '},\n'
#                    '      data: [')
#     for the_data, leg_lab in zip([large_deletions, large_insertions, coverage_array],
#                                  ['Deletions in read', 'Insertions in read', 'Total reads']):
#         html_out.write('         {\n'
#                        '         type: "line",\n'
#                        '         showInLegend: true,\n'
#                        '         legendText: "' + leg_lab + '",\n'
#                                                             '         markerType: "none",\n'
#                                                             '        dataPoints: [\n')
#         if refnum == 0:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'w')
#             bgt_out = open(
#                 options.output_folder + '/bigwig/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.bwt',
#                 'w')
#             bgt_out.write(
#                 'track type=bigwig bigDataUrl=' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                      '').lower()
#                 + '.bw name=test color=0,0,200 altColor=0,200,0 autoScale=on alwaysZero=on graphType=bar yLineMark=10 yLineOnOff=on\n')
#             bgt_out.close()
#         else:
#             bg_out = open(
#                 options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(',
#                                                                                                         '').lower() + '.wig',
#                 'a')
#         bg_out.write(
#             'fixedStep chrom=' + reference + ' start=1 step=' + str(bin_step) + ' span=' + str(bin_step) + '\n')
#         for value in range(0, len(x_axis) - 1):
#             html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
#             if x_axis[value] >= 0 and x_axis[value] < chrom_size[reference] - bin_step:
#                 bg_out.write(str(the_data[value]) + '\n')
#         bg_out.close()
#         html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
#         html_out.write('        ]\n      }')
#         if not the_data is coverage_array:
#             html_out.write(',')
#     html_out.write('''      ],
#     rangeChanged: syncHandler\n
# });
#
# chart3.render();
# var charts = [chart1, chart2, chart3];
#
# function syncHandler(e) {
#
# for (var i = 0; i < charts.length; i++) {
#     var chart = charts[i];
#
#     if (!chart.options.axisX) chart.options.axisX = {};
#
#     if (e.trigger === "reset") {
#         chart.options.axisX.viewportMinimum = chart.options.axisX.viewportMaximum = null;
#
#     } else if (chart !== e.chart) {
#         chart.options.axisX.viewportMinimum = e.axisX.viewportMinimum;
#         chart.options.axisX.viewportMaximum = e.axisX.viewportMaximum;
#     }
#
#     chart.render();
#
# }
# }
# function clickHandler(e) {
# var x = parseInt(e.target.id.split(',')[0])
# var y = parseInt(e.target.id.split(',')[1])
# for (var i = 0; i < charts.length; i++) {
#     var chart = charts[i];
#     chart.options.axisX.viewportMinimum = x;
#     chart.options.axisX.viewportMaximum = y;
#     chart.render();
#     }
# }
# var zoomButtons = document.getElementsByClassName("zoom");
# for (var i = 0; i < zoomButtons.length; i++) {
# var zoomButton = zoomButtons[i]
# zoomButton.addEventListener("click", clickHandler)
# };
# }
# </script>
# <script type="text/javascript" src="/igb/webpage_css_js/canvasjs.min.js"></script></head>
# <h1> ''' + str(html_name) + ''' graphs </h1>
# <div id="chartContainer1" style="height: 400px; width: 100%;">
# </div>
# <br />
# <table style="width: 60%" align="center">
#     <tr>
#         <th>Type</th>
#         <th>Position</th>
#     </tr>
# ''')
#     for j in clipped_flag:
#         html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(
#             j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
#         html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
#     for j in unclipped_flag:
#         html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(
#             j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
#         html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
#     if unclipped_flag == [] and clipped_flag == []:
#         html_out.write('        <tr>\n          <td> no flags </td>\n')t.write('          <td> no flags </td>\n        </tr>\n')
#     html_out.write('''
# </table>
# <br />
# <div id="chartContainer2" style="height: 400px; width: 100%;">
# </div>
# <br />
# <div id="chartContainer3" style="height: 400px; width: 100%;">
# </div>
# <br />
# <table style="width: 60%" align="center">
#     <tr>
#         <th>Type</th>
#         <th>Position</th>
#     </tr>
# ''')
#     for j in indel_flag:
#         html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(
#             j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
#         html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
#     if indel_flag == []:
#         html_out.write('        <tr>\n          <td> no flags </td>\n')
#         html_out.write('          <td> no flags </td>\n        </tr>\n')
#     html_out.write('''
# </table>
# <br />
# ''')
#     html_out.write(footer)
#     html_out.close()
#     if sam.lengths[refnum] >= 300000:
#         indexa = 100000 / bin_step
#     else:
#         indexa = sam.lengths[refnum] / 3 / bin_step
#     new_array = coverage_array[indexa:-indexa]
#     if len(new_array) >= 1000:
#         zesteps = len(new_array) / 1000
#     else:
#         zesteps = 1
#     covlist = []
#     for jump in range(0, len(new_array) + 1, zesteps):
#         covlist.append(sum(new_array[jump:jump + zesteps]) * 1.0 / zesteps)
#     out_cov[html_name] = covlist
#     out_flag[html_name] = (len(clipped_flag), len(unclipped_flag), len(indel_flag))
#     csfile = open(options.output_folder + '/chrom.size', 'w')
#     for i in chrom_size:
#         csfile.write(i + '\t' + str(chrom_size[i]) + '\n')
#     csfile.close()
#     for i in os.listdir(options.output_folder + '/wiggle/'):
#         subprocess.Popen(
#             'wigToBigWig ' + options.output_folder + '/wiggle/' + i + ' ' + options.output_folder + '/chrom.size ' + options.output_folder + '/bigwig/' + i[
#                                                                                                                                                           :-3] + 'bw',
#             shell=True).wait()
#     return out_cov, out_flag


def get_busco(busco_folder):
    bac_busco_dict = {}
    euk_busco_dict = {}
    best_busco_dict = {}
    for busco_file in os.listdir(busco_folder):
        b_path = os.path.join(busco_folder, busco_file)
        if os.path.isdir(b_path) and not '_tmp' in busco_file:
            # bin = busco_file.split('.')[1]
            bin = '.'.join(busco_file.split('.')[1:])
            summary_file_path = glob.glob(os.path.join(b_path, "short_summary*"))
            if len(summary_file_path) == 1:
                with open(summary_file_path[0]) as f:
                    for i in range(9):
                        line = f.readline()
                        if i == 8:
                            busco_string = line.split()[0]
            if busco_file.startswith('bacteria_odb10'):
                bac_busco_dict[bin] = busco_string
            elif busco_file.startswith('eukaryota_odb10'):
                euk_busco_dict[bin] = busco_string
            else:
                complete_b = float(busco_string.split(':')[1].split('%')[0])
                # TODO confirm this works? Used to cut off name in bin_summary.tsv, e.g. eria_odb10 should be bacteria_odb10
                kingdom = busco_file.split('.')[0]
                if bin in best_busco_dict:
                    if complete_b > best_busco_dict[bin][2]:
                        best_busco_dict[bin] = (kingdom, busco_string, complete_b)
                else:
                    best_busco_dict[bin] = (kingdom, busco_string, complete_b)
    return (bac_busco_dict, euk_busco_dict, best_busco_dict)


# def create_main_page(outfile, fasta, checkm_file, contig_folder, long_bam, short_bam, gff_file,
# long_qc_html, short_qc_html, gtdbtk_dir, busco_dir, instrain_file, min_contig_size=100000):
def create_main_page(outfile, fasta, checkm_file, contig_folder, long_bam, short_bam, gff_file,
                     long_qc_html, short_qc_html, gtdbtk_dir, busco_dir, min_contig_size=100000):
    with open(fasta) as f:
        len_dict = {}
        for line in f:
            if line.startswith('>'):
                name = line.split()[0][1:]
                len_dict[name] = 0
            else:
                len_dict[name] += len(line.rstrip())
    checkm_dict = {}
    with open(checkm_file) as f:
        for line in f:
            if line.startswith('-----') or line.startswith('[') or line.startswith('  Bin Id'):
                pass
            else:
                bin_id, marker_lineage, genomes, markers, marker_sets, e0, e1, e2, e3, e4, e5, completeness, contamination, heterogeneity = [
                    s.strip() for s in line.rstrip().split('  ') if s]
                checkm_dict[bin_id] = [marker_lineage, completeness, contamination, heterogeneity]
    outlist = []
    bin_list = []
    # FIXME : gene_size_dict can be missing contig entries due to incorrectly formatted or
    # missing gff entries. Either add a check below or ensure gff file contains entries for all contigs.
    gene_size_dict = get_gene_sizes(gff_file)
    instrain_dict = {}
    # with open(instrain_file) as f:
    #     f.readline()
    #     for line in f:
    #         if line.split('\t')[9] != '':
    #             instrain_dict[line.split()[0]] = float(line.split('\t')[9])
    for i in os.listdir(contig_folder):
        if not i.endswith('.fa'):
            continue
        #        bin = i.split('.')[1]
        bin_list.append(i[:-3])
    gtdbtk_dict = get_gtdbtk(gtdbtk_dir)
    busco_bac_dict, busco_euk_dict, busco_best_dict = get_busco(busco_dir)
    phylo_dict = {}
    for i in gtdbtk_dict:
        phylo_dict[i] = gtdbtk_dict[i][2]
        gtdbtk_dict[i] = gtdbtk_dict[i][0] + ' (' + gtdbtk_dict[i][1] + '%)'
    cov_dict = {}
    for i in os.listdir(contig_folder):
        if not i.endswith('.fa'):
            continue
        bin = i[:-3]  # .split('.')[1]
        ctgs = []
        with open(os.path.join(contig_folder, i)) as fa:
            for line in fa:
                if line.startswith('>'):
                    ctg_name = line.split()[0][1:]
                    ctgs.append(ctg_name)
        length_list = []
        microd = 0
        lenmicrod = 0
        for j in ctgs:
            length_list.append(len_dict[j])
            if j in instrain_dict:
                microd += len_dict[j] * instrain_dict[j]
                lenmicrod += len_dict[j]
            else:
                instrain_dict[j] = 'n/a'
        max_contig = max(length_list)
        bases_assembled = sum(length_list)
        if lenmicrod != 0:
            microd = str(microd / lenmicrod)
        else:
            microd = 'n/a'
        length_list.sort(reverse=True)
        y = 0
        for j in length_list:
            y += j
            if y >= bases_assembled / 2:
                break
        n50 = j
        ctg_details = []
        bases_sequenced_long = 0
        bases_sequenced_short = 0
        gene_sizes = []
        coding = 0
        noncoding = 0
        for ctg in ctgs:
            if ctg in gene_size_dict:
                gene_sizes += gene_size_dict[ctg]
                coding += sum(gene_size_dict[ctg])
                noncoding += len_dict[ctg] - sum(gene_size_dict[ctg])
            else:
                gene_sizes += [0]
                coding += 0
                noncoding += len_dict[ctg]
            cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in, ends_ind, x = get_cov_stats_long(
                long_bam, ctg)
            coverage_long = (sum(cov_forward) + sum(cov_reverse)) / len(cov_forward)
            bases_sequenced_long += coverage_long * len_dict[ctg]
            if short_bam is None:
                cov_forward_ill, cov_reverse_ill = None, None
                coverage_short = 0
            else:
                cov_forward_ill, cov_reverse_ill = get_cov_stats_short(short_bam, ctg)
                coverage_short = (sum(cov_forward_ill) + sum(cov_reverse_ill)) / len(cov_forward_ill)
                bases_sequenced_short += coverage_short * len_dict[ctg]
            # if len_dict[ctg] >= min_contig_size:
            #     create_contig_page(bin, ctg, cov_forward, cov_reverse, trimmed_starts, trimmed_ends, starts_in,
            #                     ends_ind, cov_forward_ill, cov_reverse_ill, x, 'ctg/' + ctg + '.html')
            ctg_details.append((ctg, len_dict[ctg], coverage_long, coverage_short, instrain_dict[ctg]))
        gene_average = numpy.average(gene_sizes)
        gene_std = numpy.std(gene_sizes)
        gene_no = len(gene_sizes)
        coding_percent = coding / (coding + noncoding) * 100 if coding != 0 else 0
        if short_bam is None:
            cov_dict[bin] = bases_sequenced_long / bases_assembled
        else:
            cov_dict[bin] = bases_sequenced_short / bases_assembled
        if bin in gtdbtk_dict:
            gtdbtk_info1 = gtdbtk_dict[bin]
            gtdbtk_info2 = phylo_dict[bin]
        else:
            gtdbtk_info1 = 'n/a'
            gtdbtk_info2 = 'n/a'
        if bin in busco_euk_dict:
            euk_busco = busco_euk_dict[bin]
        else:
            euk_busco = 'n/a'
        if bin in busco_bac_dict:
            bac_busco = busco_bac_dict[bin]
        else:
            bac_busco = 'n/a'
        if bin in busco_best_dict:
            busco_kingdom, busco_best, busco_complete = busco_best_dict[bin]
        else:
            busco_kingdom, busco_best, busco_complete = 'n/a', 'n/a', 0
        bin_headers = ['Bin', 'Max. contig (bp)', '# of contigs', 'bases assembled', 'N50', 'average read depth (long)',
                       'average read depth (short)', 'average gene size', 'Gene size Std. dev.', '# of genes',
                       'coding density (%)', 'marker lineage',
                       'completeness', 'Contamination', 'Heterozygosity', 'Microdiversity', 'Closest ref. (% ANI)',
                       'Classification', 'Bacterial busco', 'Eukaryotic busco', 'Best kingdom', 'Kingdom busco',
                       'kingdom completeness']
        hidden_headers = ['N50', 'average gene size', 'Gene size Std. dev.', '# of genes', 'coding density (%)',
                          'marker lineage', 'Classification', 'Bacterial busco', 'Eukaryotic busco', 'Best kingdom']
        bin_details = [bin, '{:,}'.format(max_contig), '{:,}'.format(len(ctgs)), '{:,}'.format(bases_assembled),
                       '{:,}'.format(n50),
                       '{:,.2f}'.format(bases_sequenced_long / bases_assembled),
                       '{:,.2f}'.format(bases_sequenced_short / bases_assembled),
                       '{:,.2f}'.format(gene_average), '{:,.2f}'.format(gene_std), '{:,}'.format(gene_no),
                       '{:,.2f}'.format(coding_percent)
                       ] + checkm_dict[bin] + [microd] + [gtdbtk_info1, gtdbtk_info2, bac_busco, euk_busco,
                                                          busco_kingdom, busco_best, '{:,.2f}'.format(busco_complete)]
        create_bin_page(bin_headers, bin_details, ctg_details, 'bin/' + bin + '.html', bin_list, long_qc_html,
                        short_qc_html)
        outlist.append(bin_details)
    get_gtdbtk(gtdbtk_dir, cov_dict)
    with open(outfile, 'w') as o:
        o.write(create_header(bin_list, '', 'index.html', long_qc_html, short_qc_html))
        o.write(add_title("Overview of assembly", "assembled into " + str(len(bin_list)) + " bins."))
        o.write(add_main("Overview of bins", "details for each bin"))
        o.write(create_table(bin_headers, outlist, hidden_headers, 'data/bin_summary.tsv'))
        o.write(end_main())
        o.write(add_footer())


def create_bin_page(headers, bin_details, ctg_details, outfile, bin_list, long_qc_html, short_qc_html):
    with open('www/' + outfile, 'w') as o:
        o.write(create_header(bin_list, '../', bin_details[0], long_qc_html, short_qc_html))
        o.write(add_title("overview of bin " + bin_details[0], "assembled into " + bin_details[2] + " contigs."))
        main_string = ''
        for i, j in zip(headers, bin_details):
            main_string += '<b>' + i + ':</b> ' + j + '<br>\n'
        ctgs = []
        for i in ctg_details:
            ctgs.append(i[0])
        o.write(add_main("Bin details:", main_string, ctgs))
        o.write(create_table(["contig", "length", "coverage long", "coverage short", "microdiversity"], ctg_details))
        o.write(end_main())
        o.write(add_footer())


checkm_file = snakemake.input.checkm_file
contig_folder = "data/das_tool_bins/das_tool_DASTool_bins/"
fasta = snakemake.input.fasta
long_html = snakemake.input.long_reads_qc_html[4:]
short_html = snakemake.input.short_reads_qc_html[4:]
gff_file = snakemake.input.genes_gff
gtdbtk_dir = snakemake.input.gtdbtk_done[:-4]
busco_dir = snakemake.input.busco_done[:-4]
# instrain_file = snakemake.input.strain_profile
try:
    os.makedirs('www/bin')
except FileExistsError:
    pass

try:
    os.makedirs('www/css')
except FileExistsError:
    pass

write_css('www/css/style.css')

if os.path.exists("data/short_reads.fastq.gz"):
    short_bam = "data/final_short.sort.bam"
    long_bam = "data/final_long.sort.bam"
else:
    long_bam = "data/final_long.sort.bam"
    short_bam = None

# create_main_page("www/index.html", fasta, checkm_file, contig_folder, long_bam, short_bam,
# gff_file, long_html, short_html, gtdbtk_dir, busco_dir, instrain_file)
create_main_page("www/index.html", fasta,
                 checkm_file, contig_folder,
                 long_bam, short_bam,
                 gff_file, long_html,
                 short_html, gtdbtk_dir, busco_dir)