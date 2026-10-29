/* Keep the primary navigation stable while Material makes room for the footer.
 *
 * Material changes the scroll wrapper's inline height as the bottom of the
 * main content enters the viewport. Without compensating its scroll position,
 * links that were visible can drop below the newly shortened wrapper and the
 * reader must move the pointer over the sidebar to reach them. Preserve the
 * wrapper's visible bottom edge by scrolling it by exactly the height lost (or
 * regained). Material remains responsible for sizing, so the sidebar still
 * cannot overlap the footer.
 */
(function () {
  "use strict";

  function mountSidebarFooterGuard() {
    const scrollwrap = document.querySelector(
      ".md-sidebar--primary .md-sidebar__scrollwrap"
    );
    const footer = document.querySelector(".md-footer");

    if (!scrollwrap || !footer || typeof ResizeObserver === "undefined") {
      return;
    }

    let previousHeight = scrollwrap.clientHeight;
    let footerWasVisible = footer.getBoundingClientRect().top < window.innerHeight;

    const observer = new ResizeObserver(function () {
      const currentHeight = scrollwrap.clientHeight;
      const heightLost = previousHeight - currentHeight;
      const footerVisible = footer.getBoundingClientRect().top < window.innerHeight;

      if ((footerVisible || footerWasVisible) && heightLost !== 0) {
        scrollwrap.scrollTop = Math.max(0, scrollwrap.scrollTop + heightLost);
      }

      previousHeight = currentHeight;
      footerWasVisible = footerVisible;
    });

    observer.observe(scrollwrap);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", mountSidebarFooterGuard, {
      once: true
    });
  } else {
    mountSidebarFooterGuard();
  }
})();
