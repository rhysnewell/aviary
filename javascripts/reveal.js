// Scroll-reveal for .aviary-reveal elements, the start signal for the
// hero's CSS-only draw-in animations, and a small typewriter cycle for the
// decorative terminal. Progressive enhancement only: every element this
// touches is fully visible/legible with this script absent (see
// extra.css's `.aviary-reveal` base state and the terminal's static
// default text), and nothing here runs under `prefers-reduced-motion: reduce`.
(function () {
  var reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  document.documentElement.classList.add('aviary-js');

  if (reduceMotion) {
    document.documentElement.classList.add('aviary-reduced-motion');
    return;
  }

  document.documentElement.classList.add('aviary-motion-ready');

  // One observer now the scrollytelling section is gone (removed 2026-08-11).
  // Elements reveal as soon as they edge into view -- a short settle as the
  // page arrives, not choreography. Anything that needs to be *noticed*
  // belongs in the layout, not in a scroll trigger.
  if ('IntersectionObserver' in window) {
    var observer = new IntersectionObserver(function (entries, self) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) {
          entry.target.classList.add('is-visible');
          self.unobserve(entry.target);
        }
      });
    }, { threshold: 0.15, rootMargin: '0px 0px -40px 0px' });

    document.querySelectorAll('.aviary-reveal').forEach(function (el) {
      observer.observe(el);
    });
  } else {
    // No observer support: just show everything, skip the choreography.
    document.querySelectorAll('.aviary-reveal').forEach(function (el) {
      el.classList.add('is-visible');
    });
  }

  /**
   * Type/erase-loop through a comma-separated `data-commands` list on an
   * element, starting from whatever static text it already contains (the
   * no-JS/reduced-motion fallback state).
   * @param {HTMLElement} el
   * @returns {void}
   */
  function typewriterLoop(el) {
    var words = (el.getAttribute('data-commands') || '')
      .split(',')
      .map(function (s) { return s.trim(); })
      .filter(Boolean);
    if (words.length < 2) return;

    var index = 0;
    var TYPE_MS = 55;
    var ERASE_MS = 30;
    var HOLD_MS = 1800;

    function type(word, done) {
      var pos = 0;
      (function step() {
        el.textContent = word.slice(0, pos);
        pos += 1;
        if (pos <= word.length) {
          setTimeout(step, TYPE_MS);
        } else {
          done();
        }
      })();
    }

    function erase(word, done) {
      var pos = word.length;
      (function step() {
        el.textContent = word.slice(0, pos);
        pos -= 1;
        if (pos >= 0) {
          setTimeout(step, ERASE_MS);
        } else {
          done();
        }
      })();
    }

    (function cycle() {
      var current = words[index % words.length];
      setTimeout(function () {
        erase(current, function () {
          index += 1;
          type(words[index % words.length], cycle);
        });
      }, HOLD_MS);
    })();
  }

  document.querySelectorAll('.aviary-terminal__type').forEach(typewriterLoop);
})();
