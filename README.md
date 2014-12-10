## Visual Computing Slides -- Curves

Part of the *Visual Computing* course given at [National University of Colombia](http://www.disi.unal.edu.co/) by [Jean Pierre Charalambos](http://otrolado.info) and [Andres Colubri](http://codeanticode.wordpress.com/).

Powered by [reveal](https://github.com/hakimel/reveal.js).

Made possible thanks to... 

<!--- a long list of students and links to their pages. To come ;) -->

## Setup

 ```sh
 $ git clone https://github.com/VisualComputing/Transformations.git
 $ git checkout gh-pages
 ```

## Folder Structure

    |-- css/
    |-- js/
    |-- plugin/
    |-- lib/
    |-- fig/
    |-- source.md
    
Refer to the [reveal folder structure](https://github.com/hakimel/reveal.js#folder-structure) for more details, and to the *Setup* below.

## Setup

External markdown and speaker notes, require that presentations run from a local web server. The following instructions will set up such a server as well as all of the development tasks needed to make edits to the slides source code.

1. Install [Node.js](http://nodejs.org/)

2. Install [Grunt](http://gruntjs.com/getting-started#installing-the-cli)

3. Navigate to the presentation folder

 ```sh
 $ cd Transformations
 ```

5. Install dependencies

 ```sh
 $ npm install
 ```

6. Edit the presentation contents using markdown in the `source.md` and adding figures as needed to the `fig/` folder.

7. Serve the presentation and monitor source files for changes

 ```sh
 $ grunt serve
 ```

8. Open <http://localhost:8000> to view your presentation

 You can change the port by using `grunt serve --port 8001`.

<!---

9. Update to upstream

 ```sh
 $ git remote add reveal.js https://github.com/hakimel/reveal.js.git
 $ git pull reveal.js master
 ```
-->
